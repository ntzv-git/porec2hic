#!/usr/bin/env bash
# ==============================================================================
#  HIFI-GUIDED PORE-C DIGESTION -> ALL-TO-ALL PSEUDO-HI-C PAIRS
# ==============================================================================
#  Usage : POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz ./porec2hic_hifi.sh
#  Tools: minimap2, seqkit, bedtools, gawk, (pigz). No Python.
#  Every variable below can be overridden from the environment.
#  Resume: each step writes its output under a temporary name and renames it
#  only on success. A step whose output file exists is skipped, so a re-run
#  after a crash resumes at the first missing output.
#
#  All intermediate files stay in Pore-C FASTQ order (bedtools map -g
#  porec.genome): no global sort is needed.
# ==============================================================================
set -Eeuo pipefail
export LC_ALL=C

# ------------------------------------------------------------------------------
# PIPELINE PARAMETERS
# ------------------------------------------------------------------------------
POREC_FQ=${POREC_FQ:?"POREC_FQ is not set"}
HIFI_FQ=${HIFI_FQ:?"HIFI_FQ is not set (HiFi FASTQ/FASTA or .mmi index)"}
OUTDIR=${OUTDIR:-porec2hic_out}
PREFIX=${PREFIX:-porec_hic}

MOTIF=${MOTIF:-CATG}              # palindromic motif, IUPAC codes allowed (GATC, GANTC, ...)
CUT_OFFSET=${CUT_OFFSET:-4}       # cut position inside the motif: NlaIII CATG^ = 4, DpnII ^GATC = 0
DUP_MOTIF=${DUP_MOTIF:-1}         # 1 = keep the motif (rebuilt by ligation) on both monomers

MIN_COV=${MIN_COV:-3}             # site PROTECTED if >= MIN_COV HiFi reads are continuous over [x-FLANK, x+FLANK]; otherwise cut
FLANK=${FLANK:-15}                # margin (bp) required on each side of the site (> motif length + alignment overshoot)
MAX_GAP=${MAX_GAP:-50}            # indel/hole (bp) from which a HiFi alignment is no longer continuous
MIN_MONO_LEN=${MIN_MONO_LEN:-$(( FLANK + ${#MOTIF} + 1 ))}  # monomers shorter than this (pseudo-monomers) are dropped before pairing
WRITE_MONOMERS=${WRITE_MONOMERS:-1}     # 1 = also write the FASTQ of kept monomers (used by the QC in README)

THREADS=${THREADS:-96}
MM2_PRESET=${MM2_PRESET:-map-ont} # lr:hq for ONT R10 Q20+ reads (minimap2 >= 2.27)
MM2_N=${MM2_N:-200}               # secondary hits = other HiFi reads of the locus (cap is PER WINDOW)
WINDOW=${WINDOW:-2000}            # Pore-C reads longer than this are aligned as overlapping windows (0 = off)
WINDOW_OVERLAP=${WINDOW_OVERLAP:-250}   # overlap (bp) between consecutive windows (> 2 x FLANK)
MM2_BATCH=${MM2_BATCH:-50G}       # HiFi index batch size (-I)
MM2_SPLIT=${MM2_SPLIT:-auto}      # --split-prefix: auto = only if the HiFi reads need several index batches
                                  # (it stores ALL alignments in temporary files until the end: huge on disk)
MM2_MAX_OCC=${MM2_MAX_OCC:-}     # ignore HiFi minimizers seen more than this many times (minimap2 -f INT).
                                  # Empty = minimap2 default (top 0.02%), which gets very high with redundant
                                  # HiFi reads / polyploids (log line "mid_occ = ..."), making step 2 slow.
                                  # Keep it well above HiFi depth x copies (e.g. 1000).
MM2_EXTRA=${MM2_EXTRA:-}

for tool in minimap2 seqkit bedtools gawk; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool not found" >&2; exit 1; }
done
if command -v pigz >/dev/null; then GZ="pigz -p 8"; else GZ="gzip -3"; fi

RC=$(echo "$MOTIF" | tr 'ACGTRYKMBDHVNacgtrykmbdhvn' 'TGCAYRMKVHDBNtgcayrmkvhdbn' | rev)
if [[ "${RC^^}" != "${MOTIF^^}" ]]; then
  echo "ERROR: motif $MOTIF is not palindromic (revcomp $RC): not supported" >&2; exit 1
fi
MOTIF_LEN=${#MOTIF}
if (( WINDOW > 0 && ( WINDOW_OVERLAP * 2 > WINDOW || WINDOW_OVERLAP <= 2 * FLANK ) )); then
  echo "ERROR: WINDOW_OVERLAP must be > 2 x FLANK and <= WINDOW / 2" >&2; exit 1
fi

# ------------------------------------------------------------------------------
# LOGGING, TIMING AND RESUME HELPERS
# ------------------------------------------------------------------------------
log() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
hms() { printf '%02d:%02d:%02d' $(( $1 / 3600 )) $(( $1 % 3600 / 60 )) $(( $1 % 60 )); }
STEP_NAME="setup"; STEP_TIMES=()
RERUN=0                                   # once a step runs, every later step runs again
done_already() { [[ $RERUN == 0 ]] && for f in "$@"; do [[ -e "$f" ]] || return 1; done; }
step_begin() {
  STEP_NAME=$1; STEP_START=$SECONDS; RERUN=1; echo; log "[$1] $2"
  case $1 in                             # cached counts derived from this step and the next ones are obsolete
    "STEP 1/5") rm -f counts.tsv ;;
    "STEP 2/5") for k in aligned_reads sites prot_cut; do kv_del counts.tsv "$k"; done ;;
    "STEP 3/5") for k in sites prot_cut; do kv_del counts.tsv "$k"; done ;;
    "STEP 4/5") kv_del counts.tsv prot_cut ;;
  esac
}
# timings.tsv keeps the duration of every step across runs (a resumed run shows
# the real total); counts.tsv caches the counts printed after each step so that
# a resumed run does not decompress large files again.
kv_get() { [[ -e "$1" ]] || return 0; awk -F'\t' -v k="$2" '$1 == k { v = $2 } END { if (v != "") print v }' "$1"; }
kv_del() { [[ -e "$1" ]] || return 0; awk -F'\t' -v k="$2" '$1 != k' "$1" > "$1.new"; mv "$1.new" "$1"; }
kv_set() { { [[ -e "$1" ]] && awk -F'\t' -v k="$2" '$1 != k' "$1"; printf '%s\t%s\n' "$2" "$3"; } > "$1.new"; mv "$1.new" "$1"; }
step_end() {
  local d=$(( SECONDS - STEP_START ))
  STEP_TIMES+=("$STEP_NAME|$d|run")
  kv_set timings.tsv "$STEP_NAME" "$d"
  log "[$STEP_NAME] done in $(hms "$d")"
}
step_skip() {
  STEP_TIMES+=("$1|$(kv_get timings.tsv "$1")|skipped")
  echo; log "[$1] $2 -> $3 already present, step skipped"
}
# count NAME COMMAND...: cached value of NAME, or run COMMAND once and cache it
count() {
  local k=$1 v; shift
  v=$(kv_get counts.tsv "$k")
  if [[ -z "$v" ]]; then v=$("$@"); kv_set counts.tsv "$k" "$v"; fi
  echo "$v"
}
trap 'log "ERROR: pipeline failed during $STEP_NAME (line $LINENO); re-run the same command to resume" >&2' ERR
kill_tree() { local c; for c in $(pgrep -P "$1"); do kill_tree "$c"; done; kill "$1" 2>/dev/null || true; }
trap 'log "INTERRUPTED during $STEP_NAME; re-run the same command to resume" >&2; trap - INT TERM ERR
      for c in $(pgrep -P $$); do kill_tree "$c"; done; exit 130' INT TERM

POREC_FQ=$(realpath "$POREC_FQ"); HIFI_FQ=$(realpath "$HIFI_FQ")
mkdir -p "$OUTDIR"
cd "$OUTDIR"

PIPELINE_START=$SECONDS
echo "======================================================================"
echo "      STARTING HIFI-GUIDED PORE-C DIGESTION PIPELINE (UNIX/BASH)"
echo "      $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
echo "Pore-C File   : $POREC_FQ"
echo "HiFi File     : $HIFI_FQ"
echo "Motif         : $MOTIF (Cut Offset: $CUT_OFFSET, motif on both monomers: $DUP_MOTIF)"
echo "Protection    : >= $MIN_COV HiFi reads continuous over site +/- ${FLANK}bp (gaps < ${MAX_GAP}bp), otherwise cut"
echo "Min monomer   : ${MIN_MONO_LEN}bp (shorter monomers removed before pairing)"
echo "Threads       : $THREADS"
echo "Output dir    : $PWD"
echo "Index batch   : $MM2_BATCH"
echo "Windows       : ${WINDOW}bp, overlap ${WINDOW_OVERLAP}bp, -N $MM2_N per window"
echo "Max occ.      : ${MM2_MAX_OCC:-minimap2 default}"
echo "======================================================================"

# ------------------------------------------------------------------------------
# STEP 1: PORE-C READ LENGTHS + OVERLAPPING WINDOWS (one pass over the FASTQ)
# ------------------------------------------------------------------------------
# porec.genome: read lengths in FASTQ order (= sort order used by bedtools).
# minimap2 caps secondary hits (-N) per QUERY. With whole reads as queries, a
# long Pore-C read (many monomers) shares one cap between all its monomers, and
# most of them would get too few HiFi reads -> false cuts. Reads longer than
# WINDOW are therefore split into windows of WINDOW bp overlapping by
# WINDOW_OVERLAP bp (the last window is aligned on the read end); shorter reads
# stay whole. Each window gets its own cap, so coverage no longer depends on
# read length. Window name = <read>_W<offset>; step 2 adds the offset back and
# merges the blocks of a HiFi read across overlapping windows.
WINDOWS=porec_windows.fa.gz
if done_already porec.genome "$WINDOWS"; then
  step_skip "STEP 1/5" "Pore-C read lengths & windows" "porec.genome $WINDOWS"
else
  step_begin "STEP 1/5" "Pore-C read lengths & ${WINDOW}bp windows (overlap ${WINDOW_OVERLAP}bp)..."
  (   # run in the background so that a stop signal is handled immediately
  trap - ERR INT TERM
  seqkit fx2tab -i -j "$THREADS" "$POREC_FQ" | \
  awk -v W="$WINDOW" -v O="$WINDOW_OVERLAP" -v GEN=porec.genome.tmp -v CNT=windows.count.tmp 'BEGIN { FS = "\t" }
    {
      L = length($2); print $1 "\t" L > GEN
      if (W <= 0 || L <= W) { print ">" $1 "_W0\n" $2; nw++; next }
      for (s = 0; s + W < L; s += W - O) { print ">" $1 "_W" s "\n" substr($2, s + 1, W); nw++ }
      print ">" $1 "_W" (L - W) "\n" substr($2, L - W + 1, W); nw++
    }
    END { print nw + 0 > CNT }' | $GZ > "$WINDOWS.tmp"
  kv_set counts.tsv windows "$(cat windows.count.tmp)"; rm -f windows.count.tmp
  mv porec.genome.tmp porec.genome
  mv "$WINDOWS.tmp" "$WINDOWS"
  ) & wait $!
  step_end
fi
read -r N_POREC TOT_BP < <(count reads_bases awk '{n++; s += $2} END {print n + 0, s + 0}' porec.genome)
log "  -> Pore-C reads input        : $N_POREC"
log "  -> Total bases / mean length : $TOT_BP bp / $(awk -v t="$TOT_BP" -v n="$N_POREC" 'BEGIN {printf "%.1f", n ? t / n : 0}') bp"
N_WIN=$(count windows bash -c "$GZ -dc '$WINDOWS' | grep -c '^>' || true")
log "  -> Alignment queries (windows) : $N_WIN"

# ------------------------------------------------------------------------------
# STEP 2: PORE-C READS (query) ALIGNED ON HIFI READS (target) -> CONTINUOUS BLOCKS
# ------------------------------------------------------------------------------
# Each Pore-C monomer aligns on the HiFi reads of its locus; -N keeps the other
# HiFi reads of the locus (secondary hits). -c gives the CIGAR: an alignment is
# split into blocks at every indel >= MAX_GAP. Window coordinates are converted
# back to read coordinates. Then, for a given HiFi read (same strand), colinear
# blocks are merged when the hole between them is < MAX_GAP on the Pore-C AND on
# the HiFi read, and their diagonals differ by < MAX_GAP (this also joins the
# overlapping blocks of two consecutive windows): a HiFi read aligned on both
# flanks of a site, with a small hole caused by sequencing errors, stays
# continuous across that site.
# Output (FASTQ order, sorted by start within each read):
#   porec_read  start  end  hifi_read
BLOCKS=hifi_blocks.bed.gz
if done_already "$BLOCKS"; then
  step_skip "STEP 2/5" "Pore-C -> HiFi alignment" "$BLOCKS"
else
  step_begin "STEP 2/5" "Aligning Pore-C reads on HiFi reads (minimap2 -c, $MM2_PRESET)..."
  (   # run in the background so that a stop signal is handled immediately
  trap - ERR INT TERM
  rm -rf mm2_tmp
  SPLIT_OPT=()
  if [[ "$MM2_SPLIT" == 1 ]]; then
    SPLIT_OPT=(--split-prefix mm2_tmp/split)
  elif [[ "$MM2_SPLIT" == auto ]]; then
    if [[ "$HIFI_FQ" == *.mmi ]]; then
      SPLIT_OPT=(--split-prefix mm2_tmp/split)       # batch count of a prebuilt index is unknown: stay safe
    else
      HIFI_BP=$(seqkit stats -T -j "$THREADS" "$HIFI_FQ" | awk 'NR == 2 { print $5 }')
      BATCH_BP=$(numfmt --from=si "${MM2_BATCH^^}")
      log "  -> HiFi bases: $HIFI_BP ; index batch: $BATCH_BP -> $(( (HIFI_BP + BATCH_BP - 1) / BATCH_BP )) batch(es)"
      (( HIFI_BP > BATCH_BP )) && SPLIT_OPT=(--split-prefix mm2_tmp/split)
    fi
  fi
  if (( ${#SPLIT_OPT[@]} )); then
    mkdir -p mm2_tmp
    log "  -> several index batches: --split-prefix (temporary alignments in mm2_tmp/, output at the end)"
  else
    log "  -> single index batch: alignments streamed directly to $BLOCKS"
  fi
  minimap2 -c -x "$MM2_PRESET" -t "$THREADS" -I "$MM2_BATCH" ${SPLIT_OPT[@]+"${SPLIT_OPT[@]}"} \
      --secondary=yes -N "$MM2_N" ${MM2_MAX_OCC:+-f "$MM2_MAX_OCC"} $MM2_EXTRA "$HIFI_FQ" "$WINDOWS" 2> minimap2.log | \
  gawk -v G="$MAX_GAP" -v LONG="[0-9]{${#MAX_GAP},}[IDN]" -v CNT=aligned.count.tmp '
    BEGIN { OFS = "\t" }
    function add(q1, q2, t1, t2) {
      if (q2 <= q1) return
      n++; HT[n] = hifi; HS[n] = st; Q1[n] = q1; Q2[n] = q2; T1[n] = t1; T2[n] = t2
    }
    function flush(   i, p, m, gq, gt) {
      if (!n) return
      # merge colinear blocks of the same HiFi read
      delete K
      for (i = 1; i <= n; i++) K[i] = sprintf("%s\t%s\t%012d", HT[i], HS[i], Q1[i])
      PROCINFO["sorted_in"] = "@val_str_asc"
      p = 0; m = 0
      for (i in K) {
        if (p && HT[i] == HT[p] && HS[i] == HS[p]) {
          gq = Q1[i] - Q2[p]
          gt = (HS[i] == "+") ? T1[i] - T2[p] : T1[p] - T2[i]
          if (gq < G && gt < G && gq - gt < G && gt - gq < G) {
            if (Q2[i] > Q2[p]) Q2[p] = Q2[i]
            if (HS[i] == "+") { if (T2[i] > T2[p]) T2[p] = T2[i] }
            else              { if (T1[i] < T1[p]) T1[p] = T1[i] }
            continue
          }
        }
        if (p) { m++; MQ1[m] = Q1[p]; MQ2[m] = Q2[p]; MT[m] = HT[p] }
        p = i
      }
      m++; MQ1[m] = Q1[p]; MQ2[m] = Q2[p]; MT[m] = HT[p]
      # output sorted by start on the Pore-C read (required by bedtools map)
      delete O
      for (i = 1; i <= m; i++) O[i] = MQ1[i]
      PROCINFO["sorted_in"] = "@val_num_asc"
      for (i in O) print cur, MQ1[i], MQ2[i], MT[i]
      n = 0; nreads++
    }
    {
      match($1, /_W[0-9]+$/)                           # window name -> read name + offset
      rd = substr($1, 1, RSTART - 1); off = substr($1, RSTART + 2) + 0
      if (rd != cur) { flush(); cur = rd }
      st = $5; hifi = $6; qs = $3 + off; qe = $4 + off
      cg = ""
      for (i = 13; i <= NF; i++) if (substr($i, 1, 5) == "cg:Z:") { cg = substr($i, 6); break }
      if (cg !~ LONG) { add(qs, qe, $8 + 0, $9 + 0); next }   # no indel >= MAX_GAP
      k = split(cg, L, /[MIDNSHP=X]/, OP)
      qp = (st == "+") ? qs : qe; tp = $8 + 0; bq = qp; bt = tp
      for (i = 1; i < k; i++) {
        len = L[i] + 0; op = OP[i]
        big = (op == "I" || op == "D" || op == "N") && len >= G
        if (big) add((st == "+") ? bq : qp, (st == "+") ? qp : bq, bt, tp)
        if (op == "M" || op == "=" || op == "X" || op == "I") qp += (st == "+") ? len : -len
        if (op == "M" || op == "=" || op == "X" || op == "D" || op == "N") tp += len
        if (big) { bq = qp; bt = tp }
      }
      add((st == "+") ? bq : qp, (st == "+") ? qp : bq, bt, tp)
    }
    END { flush(); print nreads + 0 > CNT }' | $GZ > "$BLOCKS.tmp"
  log "  -> minimap2 repeat threshold: $(grep -m1 -o 'mid_occ = [0-9]*' minimap2.log || echo 'mid_occ = ?')"
  kv_set counts.tsv aligned_reads "$(cat aligned.count.tmp)"; rm -f aligned.count.tmp
  mv "$BLOCKS.tmp" "$BLOCKS"
  rm -rf mm2_tmp
  ) & wait $!
  step_end
fi
N_ALN_READS=$(count aligned_reads bash -c "$GZ -dc '$BLOCKS' | cut -f1 | uniq | wc -l")
log "  -> Pore-C reads with HiFi alignments : $N_ALN_READS"

# ------------------------------------------------------------------------------
# STEP 3: RESTRICTION SITES (exact cut position)
# ------------------------------------------------------------------------------
# sites.tsv: read  length  motif_start  motif_end  x   (x = cut = motif_start + CUT_OFFSET)
if done_already sites.tsv || done_already sites.final.tsv.gz; then
  step_skip "STEP 3/5" "Restriction motifs" "$([[ -e sites.tsv ]] && echo sites.tsv || echo sites.final.tsv.gz)"
else
  step_begin "STEP 3/5" "Locating restriction motifs ($MOTIF)..."
  (   # run in the background so that a stop signal is handled immediately
  trap - ERR INT TERM
  seqkit locate -P -d -i -j "$THREADS" -p "$MOTIF" --bed "$POREC_FQ" | \
  awk -v off="$CUT_OFFSET" -v G=porec.genome 'BEGIN { OFS = "\t" }
    {
      while ($1 != r) {
        if ((getline line < G) <= 0) { print "ERROR: read " $1 " absent/out of order in porec.genome" > "/dev/stderr"; exit 1 }
        split(line, a, "\t"); r = a[1]; len = a[2]
      }
      x = $2 + off
      if (x > 0 && x < len) print $1, len, $2, $3, x      # x = 0 or x = length: nothing to cut
    }' > sites.tsv.tmp
  mv sites.tsv.tmp sites.tsv
  ) & wait $!
  step_end
fi
N_MOTIFS=$(count sites bash -c "if [[ -e sites.tsv ]]; then wc -l < sites.tsv; else $GZ -dc sites.final.tsv.gz | wc -l; fi")
log "  -> Candidate motif sites     : $N_MOTIFS"

# ------------------------------------------------------------------------------
# STEP 4: SITE PROTECTION BY HIFI READS (bedtools map)
# ------------------------------------------------------------------------------
#  span = number of DISTINCT HiFi reads with a continuous block covering
#         [x-FLANK, x+FLANK] ENTIRELY (window truncated at read ends)
#  P (protected): span >= MIN_COV -> not cut
#  C (cut)      : span <  MIN_COV -> cut, including when no HiFi read is present
# sites.final.tsv.gz: read  length  motif_start  motif_end  x  span  class
if done_already sites.final.tsv.gz; then
  step_skip "STEP 4/5" "HiFi protection" sites.final.tsv.gz
else
  step_begin "STEP 4/5" "Computing HiFi protection at each site (bedtools map)..."
  (   # run in the background so that a stop signal is handled immediately
  trap - ERR INT TERM
  paste sites.tsv \
    <(awk -v f="$FLANK" 'BEGIN { OFS = "\t" } { s = $5 - f; e = $5 + f; print $1, (s < 0 ? 0 : s), (e > $2 ? $2 : e) }' sites.tsv \
        | bedtools map -g porec.genome -b "$BLOCKS" -c 4 -o count_distinct -f 1.0 -a - | cut -f4) | \
  awk -v mc="$MIN_COV" 'BEGIN { OFS = "\t" }
    NF != 6 { print "ERROR: bedtools map output truncated at line " NR > "/dev/stderr"; exit 1 }
    { print $0, ($6 >= mc ? "P" : "C") }' | $GZ > sites.final.tsv.gz.tmp
  [[ $($GZ -dc sites.final.tsv.gz.tmp | wc -l) -eq $N_MOTIFS ]] || { log "ERROR: site table incomplete" >&2; exit 1; }
  mv sites.final.tsv.gz.tmp sites.final.tsv.gz
  rm -f sites.tsv
  ) & wait $!
  step_end
fi

read -r N_PROT N_CUT < <(count prot_cut bash -c "$GZ -dc sites.final.tsv.gz | awk '\$7 == \"P\" { p++ } \$7 == \"C\" { c++ } END { print p + 0, c + 0 }'")
pct() { awk -v a="$1" -v b="$2" 'BEGIN { if (b > 0) printf "%.2f", 100 * a / b; else print 0 }'; }
log "  -> Protected sites (>= $MIN_COV HiFi) : $N_PROT ($(pct "$N_PROT" "$N_MOTIFS")%)"
log "  -> Cut sites (< $MIN_COV HiFi)        : $N_CUT ($(pct "$N_CUT" "$N_MOTIFS")%)"

# ------------------------------------------------------------------------------
# STEP 5: CUTTING + PSEUDO-MONOMER FILTER + ALL-TO-ALL PAIRS
# ------------------------------------------------------------------------------
# Every unprotected site is cut. A motif lying within ~FLANK bp of a real
# junction is not protected either: a pseudo-monomer of a few bp appears between
# the two cuts. Monomers < MIN_MONO_LEN are dropped BEFORE pairing (they are not
# mappable anyway). A read with n kept monomers gives n(n-1)/2 pairs
# @read:i-j/1 and /2.
OUT5=("${PREFIX}_R1.fastq.gz" "${PREFIX}_R2.fastq.gz" "${PREFIX}.stats")
[[ "$WRITE_MONOMERS" == 1 ]] && OUT5+=("${PREFIX}_monomers.fastq.gz")
if done_already "${OUT5[@]}"; then
  step_skip "STEP 5/5" "Cutting & all-to-all pairs" "${OUT5[*]}"
else
  step_begin "STEP 5/5" "Cutting reads & writing all-to-all pseudo-Hi-C pairs..."
  (   # run in the background so that a stop signal is handled immediately
  trap - ERR INT TERM
  seqkit fx2tab -i -j "$THREADS" "$POREC_FQ" | \
  awk -v CUTS=<($GZ -dc sites.final.tsv.gz | awk '$7 == "C" { print $1 "\t" $3 "\t" $4 "\t" $5 }') \
      -v dup="$DUP_MOTIF" -v minlen="$MIN_MONO_LEN" -v wm="$WRITE_MONOMERS" \
      -v c1="$GZ > ${PREFIX}_R1.fastq.gz.tmp" -v c2="$GZ > ${PREFIX}_R2.fastq.gz.tmp" -v cm="$GZ > ${PREFIX}_monomers.fastq.gz.tmp" '
    BEGIN { FS = "\t"; more = ((getline cl < CUTS) > 0) }
    {
      id = $1; seq = $2; qual = $3; L = length(seq); nc = 0
      while (more) {
        split(cl, a, "\t")
        if (a[1] != id) break
        nc++; ncuts++; X[nc] = a[4]; MS[nc] = a[2]; ME[nc] = a[3]
        more = ((getline cl < CUTS) > 0)
      }
      n = 0; beg = 0
      for (k = 1; k <= nc + 1; k++) {
        end = (k > nc) ? L : (dup ? ME[k] : X[k])
        if (end > beg) {
          mono++
          if (end - beg < minlen) short++
          else {
            n++; sq[n] = substr(seq, beg + 1, end - beg); ql[n] = substr(qual, beg + 1, end - beg)
            if (wm) print "@" id "_" beg "_" end "\n" sq[n] "\n+\n" ql[n] | cm
          }
        }
        if (k <= nc) beg = dup ? MS[k] : X[k]
      }
      reads++; kept += n
      if (n < 2) { single++; next }
      multi++
      for (i = 1; i <= n; i++)
        for (j = i + 1; j <= n; j++) {
          tag = "@" id ":" i "-" j
          print tag "/1\n" sq[i] "\n+\n" ql[i] | c1
          print tag "/2\n" sq[j] "\n+\n" ql[j] | c2
          pairs++
        }
    }
    END {
      if (more) { print "ERROR: cut list not exhausted (" cl ") -> order mismatch with the FASTQ" > "/dev/stderr"; exit 1 }
      close(c1); close(c2); if (wm) close(cm)
      print reads + 0, mono + 0, short + 0, kept + 0, multi + 0, single + 0, pairs + 0, ncuts + 0 > "stats.tmp"
    }'
  read -r _ _ _ _ _ _ _ NCUTS_DONE < stats.tmp
  [[ $NCUTS_DONE -eq $N_CUT ]] || { log "ERROR: $NCUTS_DONE cuts applied, $N_CUT expected" >&2; exit 1; }
  for f in "${OUT5[@]}"; do [[ "$f" == "${PREFIX}.stats" ]] || mv "$f.tmp" "$f"; done
  mv stats.tmp "${PREFIX}.stats"                      # written last: marks the step as complete
  ) & wait $!
  step_end
fi

read -r READS MONO SHORT KEPT MULTI SINGLE PAIRS _ < "${PREFIX}.stats"
log "  -> Monomers generated            : $MONO"
log "  -> Monomers < ${MIN_MONO_LEN}bp removed       : $SHORT"
log "  -> Monomers kept                 : $KEPT ($(awk -v m="$KEPT" -v r="$READS" 'BEGIN {printf "%.2f", r ? m / r : 0}') monomers/read)"
log "  -> Multi-monomer reads           : $MULTI"
log "  -> Reads with < 2 monomers       : $SINGLE (no pair possible)"
log "  -> Pseudo-Hi-C pairs (all-to-all): $PAIRS"
log "  -> Output files                  : $PWD/${PREFIX}_R1.fastq.gz ($(du -h "${PREFIX}_R1.fastq.gz" | cut -f1))"
log "                                     $PWD/${PREFIX}_R2.fastq.gz ($(du -h "${PREFIX}_R2.fastq.gz" | cut -f1))"
log "  -> Per-site table                : $PWD/sites.final.tsv.gz"
log "     (read len motif_start motif_end cut span class)"
echo
log "Execution times:"
TOTAL_ALL=0; MISSING_T=0
for t in "${STEP_TIMES[@]}"; do
  IFS='|' read -r name d how <<< "$t"
  if [[ -z "$d" ]]; then log "  $name : skipped (duration not recorded)"; MISSING_T=1; continue; fi
  TOTAL_ALL=$(( TOTAL_ALL + d ))
  if [[ $how == run ]]; then log "  $name : $(hms "$d")"; else log "  $name : $(hms "$d") (previous run, skipped now)"; fi
done
log "  TOTAL (all steps) : $(hms "$TOTAL_ALL")$([[ $MISSING_T == 1 ]] && echo ' + steps without recorded duration')"
log "  TOTAL (this run)  : $(hms $(( SECONDS - PIPELINE_START )))"
echo "======================================================================"
echo "        ALL PIPELINE STEPS COMPLETED - $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
