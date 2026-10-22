#!/usr/bin/env bash
# ==============================================================================
#  FAST IN-SILICO DIGESTION OF PORE-C READS -> PSEUDO-HI-C PAIRS (ONE PASS)
# ==============================================================================
#  Usage : POREC_FQ=porec.fq.gz ./porec2hic_digest.sh
#  Tools: coreutils (split, paste), mawk (or gawk), pigz (or gzip).
#  Optional, faster decompression: rapidgzip or igzip.
#  Every variable below can be overridden from the environment.
#
#  Replaces "pore-c-py digest" (BAM) + a BAM -> FASTQ pairing script with a
#  single streaming pass, FASTQ.gz -> R1/R2.fastq.gz, with no intermediate BAM:
#    decompress | paste (1 line per read) | split -n r/W (round robin)
#      -> W mawk workers: cut at every motif, drop short monomers, write pairs
#         to their own gzip chunks
#    -> gzip chunks concatenated in the same order for R1 and R2 (a
#       concatenation of gzip files is a valid gzip file), so R1/R2 stay in sync.
#  Every motif is cut (no HiFi guidance: see porec2hic_hifi.sh for that).
# ==============================================================================
set -Eeuo pipefail
export LC_ALL=C

usage() {
  cat <<'EOF'
Usage: POREC_FQ=porec.fq.gz [VAR=value ...] porec2hic_digest.sh [-h|--help]

Cuts Pore-C reads at every restriction site and writes pseudo-Hi-C pairs
(R1/R2) in one streaming pass. Parameters are environment variables:

  POREC_FQ        Pore-C reads, FASTQ or FASTQ.gz                       (required)
  OUTDIR          output directory                                      [porec2hic_digest_out]
  PREFIX          output file prefix                                    [porec_hic]
  MOTIF           palindromic motif, IUPAC codes allowed                [CATG]
  CUT_OFFSET      cut position inside the motif (CATG^ = 4, ^GATC = 0)  [4]
  DUP_MOTIF       1 = keep the motif on both monomers                   [0]
  MIN_MONO_LEN    shorter monomers are dropped before pairing (bp)      [30]
  PAIR_MODE       all = all-to-all, n(n-1)/2 pairs; seq = adjacent only [all]
  WRITE_MONOMERS  1 = also write PREFIX_monomers.fastq.gz               [0]
  THREADS         total CPU threads                                     [96]
  WORKERS         parallel mawk workers                                 [THREADS / 4]
  GZ_THREADS      pigz threads per output stream of each worker         [2]
  GZ_LEVEL        gzip compression level (1 = fastest, 9 = smallest)    [6]
  KEEP_CHUNKS     1 = keep per-worker chunks instead of concatenating   [0]

Outputs: OUTDIR/PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz, PREFIX.stats
Pair names: @read:i-j/1 and /2 (i, j = monomer ranks in the read).

Example:
  nohup env POREC_FQ=porec.fq.gz THREADS=96 bash porec2hic_digest.sh > digest.log 2>&1 &
EOF
}
case "${1:-}" in
  -h|--help) usage; exit 0 ;;
  "") ;;
  *) echo "ERROR: unknown argument '$1' (parameters are environment variables)" >&2; usage >&2; exit 1 ;;
esac
[[ -n "${POREC_FQ:-}" ]] || { echo "ERROR: POREC_FQ is not set" >&2; usage >&2; exit 1; }

# ------------------------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------------------------
OUTDIR=${OUTDIR:-porec2hic_digest_out}
PREFIX=${PREFIX:-porec_hic}

MOTIF=${MOTIF:-CATG}              # palindromic motif, IUPAC codes allowed (GATC, GANTC, ...)
CUT_OFFSET=${CUT_OFFSET:-4}       # cut position inside the motif: NlaIII CATG^ = 4, DpnII ^GATC = 0
DUP_MOTIF=${DUP_MOTIF:-0}         # 1 = keep the motif (rebuilt by ligation) on both monomers
                                  # 0 = plain cut at CUT_OFFSET, as pore-c-py digest does
MIN_MONO_LEN=${MIN_MONO_LEN:-30}  # monomers shorter than this are dropped before pairing
PAIR_MODE=${PAIR_MODE:-all}       # all = all-to-all, n(n-1)/2 pairs ; seq = adjacent monomers only, n-1 pairs
WRITE_MONOMERS=${WRITE_MONOMERS:-0}     # 1 = also write the FASTQ of kept monomers

THREADS=${THREADS:-96}
WORKERS=${WORKERS:-$(( THREADS / 4 > 0 ? THREADS / 4 : 1 ))}  # mawk processes; each has its own compressors
GZ_THREADS=${GZ_THREADS:-2}       # pigz threads per output stream of each worker
GZ_LEVEL=${GZ_LEVEL:-6}           # 6 = pigz default; 1 is faster but gives ~30% larger files
KEEP_CHUNKS=${KEEP_CHUNKS:-0}     # 1 = do not concatenate: keep W chunks per side
                                  # (hifiasm takes them as a comma-separated list)

if command -v mawk >/dev/null; then AWK=mawk; elif command -v gawk >/dev/null; then AWK=gawk; else AWK=awk; fi
if command -v pigz >/dev/null; then GZ="pigz -$GZ_LEVEL -p $GZ_THREADS"; else GZ="gzip -$GZ_LEVEL"; fi
case "$POREC_FQ" in
  *.gz)
    if   command -v rapidgzip >/dev/null; then UNZIP="rapidgzip -dc -P 8"
    elif command -v igzip     >/dev/null; then UNZIP="igzip -dc -T 4"
    elif command -v pigz      >/dev/null; then UNZIP="pigz -dc"
    else UNZIP="gzip -dc"; fi ;;
  *) UNZIP="cat" ;;
esac

RC=$(echo "$MOTIF" | tr 'ACGTRYKMBDHVNacgtrykmbdhvn' 'TGCAYRMKVHDBNtgcayrmkvhdbn' | rev)
if [[ "${RC^^}" != "${MOTIF^^}" ]]; then
  echo "ERROR: motif $MOTIF is not palindromic (revcomp $RC): not supported" >&2; exit 1
fi
MOTIF=${MOTIF^^}
[[ "$PAIR_MODE" == all || "$PAIR_MODE" == seq ]] || { echo "ERROR: PAIR_MODE must be 'all' or 'seq'" >&2; exit 1; }
(( CUT_OFFSET >= 0 && CUT_OFFSET <= ${#MOTIF} )) || { echo "ERROR: CUT_OFFSET must be in [0, ${#MOTIF}]" >&2; exit 1; }
# IUPAC motif -> regex (a plain ACGT motif stays a plain string)
MOTIF_RE=$(echo "$MOTIF" | sed -e 's/R/[AG]/g; s/Y/[CT]/g; s/K/[GT]/g; s/M/[AC]/g; s/S/[CG]/g; s/W/[AT]/g' \
                               -e 's/B/[CGT]/g; s/D/[AGT]/g; s/H/[ACT]/g; s/V/[ACG]/g; s/N/[ACGT]/g')

log() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
hms() { printf '%02d:%02d:%02d' $(( $1 / 3600 )) $(( $1 % 3600 / 60 )) $(( $1 % 60 )); }
kill_tree() { local c; for c in $(pgrep -P "$1"); do kill_tree "$c"; done; kill "$1" 2>/dev/null || true; }
trap 'log "ERROR: digestion failed (line $LINENO)" >&2' ERR
trap 'log "INTERRUPTED" >&2; trap - INT TERM ERR; for c in $(pgrep -P $$); do kill_tree "$c"; done; exit 130' INT TERM

POREC_FQ=$(realpath "$POREC_FQ")
mkdir -p "$OUTDIR"
cd "$OUTDIR"
CHUNKS=chunks
rm -rf "$CHUNKS"; mkdir -p "$CHUNKS"

START=$SECONDS
echo "======================================================================"
echo "      PORE-C IN-SILICO DIGESTION -> PSEUDO-HI-C PAIRS (ONE PASS)"
echo "      $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
echo "Pore-C File   : $POREC_FQ"
echo "Motif         : $MOTIF (regex $MOTIF_RE, cut offset $CUT_OFFSET, motif on both monomers: $DUP_MOTIF)"
echo "Min monomer   : $MIN_MONO_LEN bp"
echo "Pairing mode  : $PAIR_MODE"
echo "Workers       : $WORKERS x ($AWK + 2 x ${GZ%% *} -p $GZ_THREADS), decompression: $UNZIP"
echo "Output Dir    : $PWD"
echo "======================================================================"

# ------------------------------------------------------------------------------
# WORKER: one input line = one read (header \t seq \t + \t qual)
# ------------------------------------------------------------------------------
# split(seq, P, motif) returns the pieces between non-overlapping motif
# occurrences, so the motif positions follow from the piece lengths, in linear
# time even for reads of hundreds of kb.
cat > worker.awk <<'AWK'
BEGIN { FS = "\t"; ml = length(motif) }
{
  id = $1; sub(/^@/, "", id); sub(/[ \t].*/, "", id)
  seq = $2; qual = $4; L = length(seq)
  np = split(seq, P, re)
  n = 0; beg = 0; pos = 0
  for (k = 1; k <= np; k++) {
    pos += length(P[k])                       # 0-based motif start (when k < np)
    if (k < np) { sites++; end = dup ? pos + ml : pos + off }
    else end = L
    if (end > beg) {
      mono++
      if (end - beg < minlen) short++
      else {
        n++; sq[n] = substr(seq, beg + 1, end - beg); ql[n] = substr(qual, beg + 1, end - beg)
        if (wm) print "@" id "_" beg "_" end "\n" sq[n] "\n+\n" ql[n] | cm
      }
    }
    if (k < np) { beg = dup ? pos : pos + off; pos += ml }
  }
  reads++; kept += n
  if (n < 2) { single++; next }
  multi++
  for (i = 1; i < n; i++)
    for (j = i + 1; j <= (all ? n : i + 1); j++) {
      tag = "@" id ":" i "-" j
      print tag "/1\n" sq[i] "\n+\n" ql[i] | c1
      print tag "/2\n" sq[j] "\n+\n" ql[j] | c2
      pairs++
    }
}
END {
  close(c1); close(c2); if (wm) close(cm)
  print reads + 0, sites + 0, mono + 0, short + 0, kept + 0, multi + 0, single + 0, pairs + 0 > stats
}
AWK

# split runs the filter once per chunk with $FILE = chunks/NNN
FILTER="$AWK -v re='$MOTIF_RE' -v motif='$MOTIF' -v off=$CUT_OFFSET -v dup=$DUP_MOTIF \
  -v minlen=$MIN_MONO_LEN -v all=$([[ $PAIR_MODE == all ]] && echo 1 || echo 0) -v wm=$WRITE_MONOMERS \
  -v c1=\"$GZ > \$FILE.R1.fastq.gz\" -v c2=\"$GZ > \$FILE.R2.fastq.gz\" -v cm=\"$GZ > \$FILE.monomers.fastq.gz\" \
  -v stats=\"\$FILE.stats\" -f worker.awk"

log "Digesting and pairing..."
( trap - ERR INT TERM
  $UNZIP "$POREC_FQ" | paste - - - - | split -n "r/$WORKERS" -d -a 4 --filter="$FILTER" - "$CHUNKS/"
) & wait $!

N_STATS=$(ls "$CHUNKS"/*.stats | wc -l)
(( N_STATS == WORKERS )) || { log "ERROR: $N_STATS worker stats for $WORKERS workers" >&2; exit 1; }
read -r READS SITES MONO SHORT KEPT MULTI SINGLE PAIRS < <(cat "$CHUNKS"/*.stats | \
  awk '{ for (i = 1; i <= NF; i++) s[i] += $i } END { for (i = 1; i <= 8; i++) printf "%d ", s[i]; print "" }')
log "Digestion done in $(hms $(( SECONDS - START )))"

# ------------------------------------------------------------------------------
# OUTPUTS: same chunk order for R1 and R2 -> pairs stay in sync
# ------------------------------------------------------------------------------
SIDES=(R1 R2); [[ "$WRITE_MONOMERS" == 1 ]] && SIDES+=(monomers)
if [[ "$KEEP_CHUNKS" == 1 ]]; then
  for s in "${SIDES[@]}"; do
    log "  -> $s chunks: $(ls "$PWD/$CHUNKS"/*."$s".fastq.gz | paste -sd,)"
  done
else
  for s in "${SIDES[@]}"; do
    out="${PREFIX}_${s}.fastq.gz"; : > "$out.tmp"
    for f in "$CHUNKS"/*."$s".fastq.gz; do cat "$f" >> "$out.tmp"; rm -f "$f"; done   # peak disk = output + 1 chunk
    mv "$out.tmp" "$out"
  done
  rm -rf "$CHUNKS"
fi
{ printf 'metric\tvalue\n'
  printf 'reads\t%s\n' "$READS"
  printf 'sites_cut\t%s\n' "$SITES"
  printf 'monomers_generated\t%s\n' "$MONO"
  printf 'monomers_short_removed\t%s\n' "$SHORT"
  printf 'monomers_kept\t%s\n' "$KEPT"
  printf 'reads_multi_monomer\t%s\n' "$MULTI"
  printf 'reads_lt2_monomers\t%s\n' "$SINGLE"
  printf 'pairs\t%s\n' "$PAIRS"
} > "${PREFIX}.stats"

log "  -> Reads                         : $READS"
log "  -> Restriction sites cut         : $SITES"
log "  -> Monomers generated            : $MONO"
log "  -> Monomers < ${MIN_MONO_LEN}bp removed       : $SHORT"
log "  -> Monomers kept                 : $KEPT ($(awk -v m="$KEPT" -v r="$READS" 'BEGIN {printf "%.2f", r ? m / r : 0}') monomers/read)"
log "  -> Multi-monomer reads           : $MULTI"
log "  -> Reads with < 2 monomers       : $SINGLE (no pair possible)"
log "  -> Pseudo-Hi-C pairs ($PAIR_MODE)       : $PAIRS"
if [[ "$KEEP_CHUNKS" != 1 ]]; then
  log "  -> Output files                  : $PWD/${PREFIX}_R1.fastq.gz ($(du -h "${PREFIX}_R1.fastq.gz" | cut -f1))"
  log "                                     $PWD/${PREFIX}_R2.fastq.gz ($(du -h "${PREFIX}_R2.fastq.gz" | cut -f1))"
fi
log "  TOTAL : $(hms $(( SECONDS - START )))"
echo "======================================================================"
