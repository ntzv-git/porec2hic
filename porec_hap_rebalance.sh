#!/usr/bin/env bash
# ==============================================================================
#  PORE-C-GUIDED HAPLOTYPE REBALANCING OF A PHASED POLYPLOID ASSEMBLY
# ==============================================================================
#  Usage : HAPS="hap1.fa hap2.fa ..." MONOMERS=porec_hic_monomers.fastq.gz \
#          ORG_FA=organelles.fa ./porec_hap_rebalance.sh
#  Tools: minimap2, compleasm, gawk or mawk, (pigz).
#
#  A phased assembler can put contigs of one haplotype into another one: the
#  receiving haplotype gets duplicated genes, the other one misses them.
#  Purging the duplicates would then delete the missing haplotype's sequence.
#  This script moves such contigs back, using Pore-C contacts:
#    1. organelle filter
#    2. check (compleasm): sizes, S/D/F/M, and whether the genes MISSING in one
#       haplotype are DUPLICATED in another one more often than by chance
#       -> verdict; the script stops here unless the verdict is MISASSIGNED
#    3. Pore-C monomers aligned ONCE on all haplotypes (minimap2, MAPQ filter)
#    4. per contig: number of Pore-C reads linking it to other contigs of each
#       haplotype -> a contig mostly linked to another haplotype is moved there
#    5. check again -> verdict for the next step (scaffolding or purge)
#  Resume: a step whose output exists is skipped (delete it to redo the step).
# ==============================================================================
set -Eeuo pipefail
export LC_ALL=C

usage() {
  cat <<'EOF'
Usage: HAPS="hap1.fa hap2.fa ..." MONOMERS=monomers.fastq.gz [VAR=value ...] porec_hap_rebalance.sh [-h|--help]

Moves contigs that a phased assembler put in the wrong haplotype back to the
haplotype they share Pore-C contacts with. Checks with compleasm before and
after, and writes a verdict that tells what to do next.

  HAPS            haplotype FASTA files, in haplotype order (hap1 hap2 ...)     (required)
  MONOMERS        Pore-C monomers, FASTQ(.gz), from                             (required)
                  porec2hic_digest.sh PAIR_MODE=none (names read_start_end)
  ORG_FA          organelle sequences (mitochondrion + chloroplast) FASTA       [none: no filter]
  OUTDIR          output directory                                              [hap_rebalance_out]
  THREADS         CPU threads                                                   [32]
  HAP_SIZE_MB     expected size of one haplotype (Mb), only for the size check  [mean of the haplotypes]
  LINEAGE         compleasm lineage                                             [liliopsida]
  COMPLEASM_LIB   compleasm lineage directory (-L)                              [compleasm default]

  Organelle filter: a contig is removed if all of these hold
  MIN_ORG_ALN     aligned bases on organelles                                   [30000]
  MIN_ORG_COV     fraction of the contig aligned on organelles                  [0.80]
  MIN_ORG_ID      identity of these alignments                                  [0.80]
  MAX_ORG_LEN     contigs longer than this are never removed                    [300000]

  Contacts and moves
  MM2_PRESET      minimap2 preset for the monomers (lr:hq for ONT Q20+)         [map-ont]
  MIN_MAPQ        monomers with a lower MAPQ are ignored                        [20]
  MIN_CTG_LEN     shorter contigs are never moved                               [50000]
  MIN_LINKS       Pore-C reads needed between a contig and the target haplotype [20]
  MIN_SHARE       fraction of the contig's links that go to the target          [0.5]
  MIN_RATIO       target links / links to its current haplotype                 [2]
  MAX_ROUNDS      rounds of moves (round 1: contigs carrying the flagged genes) [3]

  Verdict
  D_EXCESS        duplicated genes (points of %) above the lowest haplotype     [5]
                  that count as an excess
  ENRICH          (missing in A and duplicated in B) / expected by chance       [1.3]
                  that counts as a misassignment from A to B
  MIN_OBS         minimum genes behind such an enrichment                       [20]
  FORCE_REBALANCE 1 = run steps 3-5 whatever the first verdict                  [0]

Outputs (OUTDIR): final/hapN.fa (prefix hapN_ on every contig), DECISIONS.txt,
  2_check_before/ and 5_check_after/ (diagnosis.tsv, enrichment.tsv, verdict.txt),
  4_moves/moves.tsv and contig_links.tsv, 1_organelles/organelle_contigs.tsv
EOF
}
case "${1:-}" in
  -h|--help) usage; exit 0 ;;
  "") ;;
  *) echo "ERROR: unknown argument '$1' (parameters are environment variables)" >&2; usage >&2; exit 1 ;;
esac
[[ -n "${HAPS:-}" ]] || { echo "ERROR: HAPS is not set" >&2; usage >&2; exit 1; }
[[ -n "${MONOMERS:-}" ]] || { echo "ERROR: MONOMERS is not set" >&2; usage >&2; exit 1; }

# ------------------------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------------------------
ORG_FA=${ORG_FA:-}
OUTDIR=${OUTDIR:-hap_rebalance_out}
THREADS=${THREADS:-32}
HAP_SIZE_MB=${HAP_SIZE_MB:-}
LINEAGE=${LINEAGE:-liliopsida}
COMPLEASM_LIB=${COMPLEASM_LIB:-}
MIN_ORG_ALN=${MIN_ORG_ALN:-30000}
MIN_ORG_COV=${MIN_ORG_COV:-0.80}
MIN_ORG_ID=${MIN_ORG_ID:-0.80}
MAX_ORG_LEN=${MAX_ORG_LEN:-300000}
MM2_PRESET=${MM2_PRESET:-map-ont}
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_CTG_LEN=${MIN_CTG_LEN:-50000}
MIN_LINKS=${MIN_LINKS:-20}
MIN_SHARE=${MIN_SHARE:-0.5}
MIN_RATIO=${MIN_RATIO:-2}
MAX_ROUNDS=${MAX_ROUNDS:-3}
D_EXCESS=${D_EXCESS:-5}
ENRICH=${ENRICH:-1.3}
MIN_OBS=${MIN_OBS:-20}
FORCE_REBALANCE=${FORCE_REBALANCE:-0}

if command -v gawk >/dev/null; then AWK=gawk; elif command -v mawk >/dev/null; then AWK=mawk; else AWK=awk; fi
if command -v pigz >/dev/null; then GZ="pigz -p 8"; else GZ="gzip"; fi

log() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
hms() { printf '%02d:%02d:%02d' $(( $1 / 3600 )) $(( $1 % 3600 / 60 )) $(( $1 % 60 )); }
kill_tree() { local c; for c in $(pgrep -P "$1"); do kill_tree "$c"; done; kill "$1" 2>/dev/null || true; }
STEP_NAME="setup"
trap 'log "ERROR: failed during $STEP_NAME (line $LINENO); re-run the same command to resume" >&2' ERR
trap 'log "INTERRUPTED during $STEP_NAME; re-run the same command to resume" >&2; trap - INT TERM ERR
      for c in $(pgrep -P $$); do kill_tree "$c"; done; exit 130' INT TERM
step() { STEP_NAME=$1; echo; log "[$1] $2"; }

read -r -a HAP_FILES <<< "$HAPS"
NHAP=${#HAP_FILES[@]}
(( NHAP >= 2 )) || { echo "ERROR: HAPS must list at least 2 haplotype FASTA files" >&2; exit 1; }
for i in "${!HAP_FILES[@]}"; do
  [[ -s "${HAP_FILES[$i]}" ]] || { echo "ERROR: ${HAP_FILES[$i]} not found" >&2; exit 1; }
  HAP_FILES[$i]=$(realpath "${HAP_FILES[$i]}")
done
MONOMERS=$(realpath "$MONOMERS")
[[ -z "$ORG_FA" ]] || ORG_FA=$(realpath "$ORG_FA")
for tool in minimap2 compleasm; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool not found" >&2; exit 1; }
done

mkdir -p "$OUTDIR"; cd "$OUTDIR"
START=$SECONDS
echo "======================================================================"
echo "      PORE-C-GUIDED HAPLOTYPE REBALANCING"
echo "      $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
for i in "${!HAP_FILES[@]}"; do echo "hap$(( i + 1 ))          : ${HAP_FILES[$i]}"; done
echo "Monomers      : $MONOMERS"
echo "Organelles    : ${ORG_FA:-none (no filter)}"
echo "Moves         : contig >= $MIN_CTG_LEN bp, >= $MIN_LINKS reads, share >= $MIN_SHARE, ratio >= $MIN_RATIO, MAPQ >= $MIN_MAPQ"
echo "Output Dir    : $PWD"
echo "======================================================================"

# ------------------------------------------------------------------------------
# STEP 1: haplotype prefixes + organelle filter
# ------------------------------------------------------------------------------
# Every contig is renamed hapN_<name> (unless it already has this prefix), so
# that the haplotype of a contig can be read from its name everywhere.
step "STEP 1/5" "Haplotype prefixes and organelle filter"
mkdir -p 1_organelles 1_haps
for i in "${!HAP_FILES[@]}"; do
  h=hap$(( i + 1 ))
  [[ -s 1_haps/$h.fa ]] && continue
  $AWK -v p="${h}_" '/^>/ { n = substr($1, 2); if (index(n, p) != 1) n = p n; print ">" n; next } { print }' \
    "${HAP_FILES[$i]}" > 1_haps/$h.fa.tmp
  mv 1_haps/$h.fa.tmp 1_haps/$h.fa
done
if [[ ! -s 1_organelles/organelle_contigs.tsv ]]; then
  printf 'contig\tlength\torg_aligned\torg_cov\torg_id\n' > 1_organelles/organelle_contigs.tsv.tmp
  if [[ -n "$ORG_FA" ]]; then
    cat 1_haps/hap*.fa > 1_organelles/all_haps.fa
    minimap2 -t "$THREADS" -x asm20 "$ORG_FA" 1_organelles/all_haps.fa 2> 1_organelles/minimap2.log | \
    sort -k1,1 -k3,3n | $AWK -F'\t' -v OFS='\t' -v mal="$MIN_ORG_ALN" -v mcov="$MIN_ORG_COV" -v mid="$MIN_ORG_ID" -v mlen="$MAX_ORG_LEN" '
      function flush() {
        if (q == "") return
        al = u + (e - s); cov = al / ql; id = alen ? m / alen : 0
        if (al >= mal && cov >= mcov && id >= mid && ql <= mlen) printf "%s\t%d\t%d\t%.3f\t%.3f\n", q, ql, al, cov, id
      }
      $1 != q { flush(); q = $1; ql = $2; u = 0; s = $3; e = $4; m = 0; alen = 0 }
      { m += $10; alen += $11
        if ($3 > e) { u += e - s; s = $3; e = $4 } else if ($4 > e) e = $4 }   # union of query intervals
      END { flush() }' >> 1_organelles/organelle_contigs.tsv.tmp
    rm -f 1_organelles/all_haps.fa
  fi
  mv 1_organelles/organelle_contigs.tsv.tmp 1_organelles/organelle_contigs.tsv
fi
mkdir -p 1_haps_noorg
for i in $(seq 1 "$NHAP"); do
  [[ -s 1_haps_noorg/hap$i.fa ]] && continue
  $AWK 'NR == FNR { if (FNR > 1) drop[$1]; next }
        /^>/ { keep = !(substr($1, 2) in drop) } keep' 1_organelles/organelle_contigs.tsv 1_haps/hap$i.fa > 1_haps_noorg/hap$i.fa.tmp
  mv 1_haps_noorg/hap$i.fa.tmp 1_haps_noorg/hap$i.fa
done
log "  -> organelle contigs removed: $(( $(wc -l < 1_organelles/organelle_contigs.tsv) - 1 )) (1_organelles/organelle_contigs.tsv)"

# ------------------------------------------------------------------------------
# CHECK: compleasm + sizes + enrichment -> verdict
# ------------------------------------------------------------------------------
# diagnose DIR FASTA_DIR: writes DIR/diagnosis.tsv, DIR/enrichment.tsv,
# DIR/genes.tsv (hap, gene, status, contig) and DIR/verdict.txt
diagnose() {
  local dir=$1 fadir=$2 i
  mkdir -p "$dir/compleasm"
  for i in $(seq 1 "$NHAP"); do
    if ! find "$dir/compleasm/hap$i" -name full_table_busco_format.tsv 2>/dev/null | grep -q .; then
      log "  compleasm on hap$i..."
      compleasm run -a "$fadir/hap$i.fa" -o "$dir/compleasm/hap$i" -l "$LINEAGE" -t "$THREADS" \
        ${COMPLEASM_LIB:+-L "$COMPLEASM_LIB"} > "$dir/compleasm/hap$i.log" 2>&1
    fi
  done
  for i in $(seq 1 "$NHAP"); do
    $AWK -F'\t' -v OFS='\t' -v h="hap$i" '!/^#/ && $1 != "" { print h, $1, $2, $3 }' \
      "$(find "$dir/compleasm/hap$i" -name full_table_busco_format.tsv | head -n 1)"
  done > "$dir/genes.tsv"
  for i in $(seq 1 "$NHAP"); do
    $AWK -v h="hap$i" '/^>/ { n++; next } { t += length($0) } END { printf "%s\t%d\t%.1f\n", h, n, t / 1e6 }' "$fadir/hap$i.fa"
  done > "$dir/sizes.tsv"

  # gene status per haplotype: S (Single = "Complete"), D, F, M
  $AWK -F'\t' -v OFS='\t' -v hs="$HAP_SIZE_MB" '
    NR == FNR { n[$1] = $2; mb[$1] = $3; tot += $3; nh++; next }
    { k = $1 SUBSEP $2; if (!(k in st) || $3 == "Duplicated") st[k] = $3; g[$2] }
    END {
      ng = length(g); exp_mb = (hs != "") ? hs : tot / nh
      for (k in st) { split(k, a, SUBSEP); c[a[1], st[k]]++ }
      print "hap", "contigs", "size_Mb", "size_ratio", "S%", "D%", "F%", "M%"
      for (i = 1; i <= nh; i++) {
        h = "hap" i
        printf "%s\t%d\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", h, n[h], mb[h], mb[h] / exp_mb,
          100 * c[h, "Complete"] / ng, 100 * c[h, "Duplicated"] / ng, 100 * c[h, "Fragmented"] / ng, 100 * c[h, "Missing"] / ng
      }
    }' "$dir/sizes.tsv" "$dir/genes.tsv" > "$dir/diagnosis.tsv"

  # enrichment: genes missing in A and duplicated in B, observed vs expected by chance
  $AWK -F'\t' -v OFS='\t' -v enr="$ENRICH" -v mobs="$MIN_OBS" '
    { k = $1 SUBSEP $2; if (!(k in st) || $3 == "Duplicated") st[k] = $3; g[$2]; h[$1] }
    END {
      ng = length(g)
      for (k in st) { split(k, a, SUBSEP); if (st[k] == "Missing") { M[a[1], a[2]]; nm[a[1]]++ } if (st[k] == "Duplicated") { D[a[1], a[2]]; nd[a[1]]++ } }
      print "missing_in", "duplicated_in", "observed", "expected", "ratio", "flag"
      for (A in h) for (B in h) if (A != B && nm[A] > 0) {
        o = 0; for (x in g) if (((A, x) in M) && ((B, x) in D)) o++
        e = nm[A] * nd[B] / ng; r = e > 0 ? o / e : 0
        printf "%s\t%s\t%d\t%.1f\t%.2f\t%s\n", A, B, o, e, r, (r >= enr && o >= mobs && o - e >= mobs / 2) ? "MISASSIGNED" : "-"
      }
    }' "$dir/genes.tsv" | { read -r hdr; echo "$hdr"; sort -t$'\t' -k5,5gr; } > "$dir/enrichment.tsv"

  # verdict
  $AWK -F'\t' -v dex="$D_EXCESS" '
    NR == FNR { if (FNR > 1) { d[$1] = $6; m[$1] = $8; if (!nd++ || $6 < dmin) dmin = $6 } next }
    FNR > 1 && $6 == "MISASSIGNED" { mis = mis "\n  " $2 " holds contigs of " $1 " (" $3 " genes missing in " $1 " are duplicated in " $2 ", " $4 " expected)" }
    END {
      for (h in d) if (d[h] > dmin + dex) red = red " " h
      for (h in m) if (m[h] > 5) inc = inc " " h
      if (mis != "") {
        print "VERDICT: MISASSIGNED" mis
        print "NEXT: move these contigs back with Pore-C contacts (steps 3-5). Do NOT purge duplicates before that."
        if (red != "") print "NOTE: duplicated genes also above the lowest haplotype + " dex " points in:" red " (see the check after rebalancing)"
      } else if (red != "") {
        print "VERDICT: REDUNDANT (duplicated genes above the lowest haplotype + " dex " points in:" red ")"
        print "NEXT: no misassignment detected: run purge_dups on these haplotypes only, then scaffold."
      } else {
        print "VERDICT: BALANCED"
        print "NEXT: no purge needed: go to scaffolding."
      }
      if (inc != "") print "WARNING: more than 5% missing genes in:" inc " (not explained by another haplotype if no MISASSIGNED line names them)"
    }' "$dir/diagnosis.tsv" "$dir/enrichment.tsv" > "$dir/verdict.txt"
}

show_check() {
  local dir=$1
  log "  diagnosis ($dir/diagnosis.tsv):"; sed 's/^/      /' "$dir/diagnosis.tsv"
  log "  strongest enrichments ($dir/enrichment.tsv):"; head -n 6 "$dir/enrichment.tsv" | sed 's/^/      /'
  log "  $(head -n 1 "$dir/verdict.txt")"; tail -n +2 "$dir/verdict.txt" | sed 's/^/      /'
}

step "STEP 2/5" "Check before rebalancing (compleasm $LINEAGE)"
[[ -s 2_check_before/verdict.txt ]] || diagnose 2_check_before 1_haps_noorg
show_check 2_check_before
{ echo "== before rebalancing"; cat 2_check_before/verdict.txt; } > DECISIONS.txt

if ! grep -q "VERDICT: MISASSIGNED" 2_check_before/verdict.txt && [[ "$FORCE_REBALANCE" != 1 ]]; then
  mkdir -p final; cp 1_haps_noorg/hap*.fa final/
  log "No misassignment: steps 3-5 skipped, final/ = haplotypes without organelles (FORCE_REBALANCE=1 to run them anyway)"
  log "TOTAL : $(hms $(( SECONDS - START )))"
  exit 0
fi

# ------------------------------------------------------------------------------
# STEP 3: monomers aligned once on all haplotypes
# ------------------------------------------------------------------------------
# minimap2 keeps the input order, so the monomers of a read stay together.
# One line per kept monomer: read <tab> contig (best alignment, MAPQ >= MIN_MAPQ).
step "STEP 3/5" "Aligning Pore-C monomers on all haplotypes (minimap2 -x $MM2_PRESET)"
mkdir -p 3_contacts
if [[ ! -s 3_contacts/monomer_hits.tsv.gz ]]; then
  cat 1_haps_noorg/hap*.fa > 3_contacts/all_haps.fa
  minimap2 -t "$THREADS" -x "$MM2_PRESET" --secondary=no 3_contacts/all_haps.fa "$MONOMERS" 2> 3_contacts/minimap2.log | \
  $AWK -F'\t' -v OFS='\t' -v mq="$MIN_MAPQ" '
    function flush() { if (q != "" && bq >= mq) { r = q; sub(/_[0-9]+_[0-9]+$/, "", r); print r, bt } }
    $1 != q { flush(); q = $1; bm = -1 }
    $10 > bm { bm = $10; bt = $6; bq = $12 }
    END { flush() }' | $GZ > 3_contacts/monomer_hits.tsv.gz.tmp
  mv 3_contacts/monomer_hits.tsv.gz.tmp 3_contacts/monomer_hits.tsv.gz
  rm -f 3_contacts/all_haps.fa
fi
log "  -> monomers placed (MAPQ >= $MIN_MAPQ): $($GZ -dc 3_contacts/monomer_hits.tsv.gz | wc -l)"

# ------------------------------------------------------------------------------
# STEP 4: contig links, moves
# ------------------------------------------------------------------------------
# pairs(C, D) = number of Pore-C reads with monomers on both contigs C and D.
# links(C, H) = sum of pairs(C, D) over the contigs D currently in haplotype H.
# Moves are only allowed from a haplotype B to a haplotype A that the check
# flagged as "A misses genes duplicated in B" (MISASSIGNED).
#   round 1: candidates = contigs of B carrying such genes; their links are
#            counted with non-candidate contigs only (anchors), because the
#            misplaced contigs of a region are often misplaced together.
#   round 2+: every other contig of B, links with the corrected assignment;
#            catches misplaced contigs without compleasm genes. Stops when a
#            round moves nothing (at most MAX_ROUNDS rounds).
# A contig moves from B to A if: length >= MIN_CTG_LEN, links(C, A) >= MIN_LINKS,
# links(C, A) >= MIN_RATIO x links(C, B), and links(C, A) >= MIN_SHARE x all its links.
step "STEP 4/5" "Counting Pore-C links between contigs, choosing moves"
mkdir -p 4_moves
if [[ ! -s 4_moves/contig_pairs.tsv.gz ]]; then
  $GZ -dc 3_contacts/monomer_hits.tsv.gz | $AWK -F'\t' -v OFS='\t' '
    function flush(   i, j) {
      for (i = 1; i < n; i++) for (j = i + 1; j <= n; j++) P[C[i] < C[j] ? C[i] SUBSEP C[j] : C[j] SUBSEP C[i]]++
      delete has; n = 0
    }
    $1 != r { flush(); r = $1 }
    !($2 in has) { has[$2]; C[++n] = $2 }
    END { flush(); for (k in P) { split(k, a, SUBSEP); print a[1], a[2], P[k] } }' | $GZ > 4_moves/contig_pairs.tsv.gz.tmp
  mv 4_moves/contig_pairs.tsv.gz.tmp 4_moves/contig_pairs.tsv.gz
fi
if [[ ! -s 4_moves/moves.tsv ]]; then
  for i in $(seq 1 "$NHAP"); do
    $AWK -v OFS='\t' '/^>/ { if (n != "") print n, l; n = substr($1, 2); l = 0; next } { l += length($0) } END { if (n != "") print n, l }' 1_haps_noorg/hap$i.fa
  done > 4_moves/contig_lengths.tsv
  $GZ -dc 4_moves/contig_pairs.tsv.gz | $AWK -F'\t' -v OFS='\t' -v nh="$NHAP" -v mlen="$MIN_CTG_LEN" -v ml="$MIN_LINKS" \
      -v msh="$MIN_SHARE" -v mr="$MIN_RATIO" -v maxr="$MAX_ROUNDS" -v links_out=4_moves/contig_links.tsv '
    function hapof(c) { return substr(c, 1, index(c, "_") - 1) }
    FILENAME == ARGV[1] { len[$1] = $2; cur[$1] = hapof($1); next }
    FILENAME == ARGV[2] { if (FNR > 1 && $6 == "MISASSIGNED") { tgt[$2, $1]; donor[$2] } next }
    FILENAME == ARGV[3] { k = $1 SUBSEP $2; if (!(k in st) || $3 == "Duplicated") st[k] = $3
                          if ($3 == "Duplicated" && $4 != "") dupon[$1, $2, $4]; next }
    { nb[$1] = nb[$1] " " $2; nb[$2] = nb[$2] " " $1; w[$1, $2] = $3; w[$2, $1] = $3 }
    function links(c, anchors_only,   n, i, d, a) {
      delete L; n = split(nb[c], a, " ")
      for (i = 1; i <= n; i++) { d = a[i]; if (anchors_only && (d in cand)) continue; L[cur[d]] += w[c, d] }
    }
    END {
      for (k in dupon) { split(k, a, SUBSEP)
        for (t = 1; t <= nh; t++) if (((a[1], "hap" t) in tgt) && st["hap" t, a[2]] == "Missing") { cand[a[3]]; sup[a[3], "hap" t]++ } }
      for (r = 1; r <= maxr; r++) {
        nm = 0
        for (c in len) {
          B = cur[c]
          if (!(B in donor) || (c in moved) || len[c] < mlen || (r == 1 && !(c in cand))) continue
          links(c, r == 1)
          tot = 0; for (h in L) tot += L[h]
          A = ""; for (t = 1; t <= nh; t++) if (((B, "hap" t) in tgt) && (A == "" || L["hap" t] > L[A])) A = "hap" t
          if (A == "" || L[A] < ml || L[A] < mr * L[B] || L[A] < msh * tot) continue
          to[c] = A; nm++
          printf "%s\t%s\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\n", c, B, A, len[c], L[B], L[A], L[A] / tot, sup[c, A], r
        }
        for (c in to) { cur[c] = to[c]; moved[c] }
        delete to
        if (!nm) break
      }
      printf "contig\thap_before\thap_after\tlength\ttotal" > links_out; for (t = 1; t <= nh; t++) printf "\tlinks_hap%d", t > links_out; print "" > links_out
      for (c in len) { links(c, 0); tot = 0; for (h in L) tot += L[h]
        printf "%s\t%s\t%s\t%d\t%d", c, hapof(c), cur[c], len[c], tot > links_out
        for (t = 1; t <= nh; t++) printf "\t%d", L["hap" t] > links_out; print "" > links_out }
    }' 4_moves/contig_lengths.tsv 2_check_before/enrichment.tsv 2_check_before/genes.tsv - | sort -t$'\t' -k9,9n -k2,2 -k3,3 -k4,4nr | \
    { printf 'contig\tfrom\tto\tlength\tlinks_from\tlinks_to\tshare_to\tgenes_dup_here_missing_there\tround\n'; cat; } > 4_moves/moves.tsv.tmp
  mv 4_moves/moves.tsv.tmp 4_moves/moves.tsv
fi
log "  -> contigs moved: $(( $(wc -l < 4_moves/moves.tsv) - 1 )) (4_moves/moves.tsv; links per contig: 4_moves/contig_links.tsv)"
$AWK -F'\t' 'FNR > 1 { k = $2 " -> " $3 " (round " $9 ")"; n[k]++; bp[k] += $4; g[k] += $8 }
  END { for (k in n) printf "      %s : %d contigs, %.1f Mb, %d supporting genes\n", k, n[k], bp[k] / 1e6, g[k] }' 4_moves/moves.tsv | sort

mkdir -p 4_moves/haps
if [[ ! -s 4_moves/haps/hap$NHAP.fa ]]; then
  rm -f 4_moves/haps/*.fa
  cat 1_haps_noorg/hap*.fa | $AWK -v d=4_moves/haps 'NR == FNR { if (FNR > 1) to[$1] = $3; next }
    /^>/ { n = substr($1, 2); h = substr(n, 1, index(n, "_") - 1); if (n in to) h = to[n]
           out = d "/" h ".fa"; print ">" h substr(n, index(n, "_")) > out; next }
    { print > out }' 4_moves/moves.tsv -
fi

# ------------------------------------------------------------------------------
# STEP 5: check after rebalancing
# ------------------------------------------------------------------------------
step "STEP 5/5" "Check after rebalancing (compleasm $LINEAGE)"
[[ -s 5_check_after/verdict.txt ]] || diagnose 5_check_after 4_moves/haps
show_check 5_check_after
log "  change (after - before), points of %:"
paste 2_check_before/diagnosis.tsv 5_check_after/diagnosis.tsv | $AWK -F'\t' 'NR == 1 { print "      hap\tsize_Mb\tS%\tD%\tM%"; next }
  { printf "      %s\t%+.1f\t%+.2f\t%+.2f\t%+.2f\n", $1, $11 - $3, $13 - $5, $14 - $6, $16 - $8 }'
{ echo; echo "== after rebalancing ($(( $(wc -l < 4_moves/moves.tsv) - 1 )) contigs moved)"; cat 5_check_after/verdict.txt; } >> DECISIONS.txt
mkdir -p final; cp 4_moves/haps/hap*.fa final/
log "Final haplotypes: $PWD/final/hap*.fa ; decisions: $PWD/DECISIONS.txt"
log "TOTAL : $(hms $(( SECONDS - START )))"
echo "======================================================================"
