#!/usr/bin/env bash
# ==============================================================================
#  HIFI-GUIDED PORE-C DIGESTION -> ALL-TO-ALL PSEUDO-HI-C PAIRS
# ==============================================================================
#  Usage : POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz ./porec2hic_hifi.sh
#  Outils : minimap2, seqkit, bedtools, gawk, (pigz). Aucun script python.
#  Toutes les variables ci-dessous peuvent être surchargées par l'environnement.
#  Chaque étape laisse un fichier stepN.done : une relance reprend là où elle s'est arrêtée.
#
#  Tous les fichiers intermédiaires restent dans l'ordre du FASTQ Pore-C
#  (bedtools map -g porec.genome) : aucun tri global n'est nécessaire.
# ==============================================================================
set -euo pipefail
export LC_ALL=C

# ------------------------------------------------------------------------------
# VARIABLES DU PIPELINE
# ------------------------------------------------------------------------------
POREC_FQ=${POREC_FQ:?"POREC_FQ non défini"}
HIFI_FQ=${HIFI_FQ:?"HIFI_FQ non défini (FASTQ/FASTA HiFi ou index .mmi)"}
OUTDIR=${OUTDIR:-porec2hic_out}
PREFIX=${PREFIX:-porec_hic}

MOTIF=${MOTIF:-CATG}              # motif palindromique, IUPAC accepté (GATC, GANTC, ...)
CUT_OFFSET=${CUT_OFFSET:-4}       # coupure dans le motif : NlaIII CATG^ = 4, DpnII ^GATC = 0
DUP_MOTIF=${DUP_MOTIF:-1}         # 1 = le motif reconstitué à la ligation est gardé sur les 2 monomères

MIN_COV=${MIN_COV:-3}             # site PROTÉGÉ si >= MIN_COV reads HiFi continus sur [x-FLANK, x+FLANK] ; sinon coupé
FLANK=${FLANK:-15}                # marge (pb) exigée de part et d'autre du site (> longueur motif + débordement)
MAX_GAP=${MAX_GAP:-50}            # indel/trou (pb) au-delà duquel un alignement HiFi n'est plus continu
MIN_MONO_LEN=${MIN_MONO_LEN:-$(( FLANK + ${#MOTIF} + 1 ))}  # < FLANK+motif : pseudo-monomères écartés avant l'appariement
WRITE_MONOMERS=${WRITE_MONOMERS:-0}     # 1 = écrire aussi le FASTQ des monomères

THREADS=${THREADS:-96}
MM2_PRESET=${MM2_PRESET:-map-ont} # lr:hq pour ONT R10 Q20+ (minimap2 >= 2.27)
MM2_N=${MM2_N:-100}               # secondaires = autres reads HiFi du locus (plafond PAR READ Pore-C)
MM2_BATCH=${MM2_BATCH:-50G}       # lots d'index HiFi (-I) + --split-prefix
MM2_EXTRA=${MM2_EXTRA:-}

for tool in minimap2 seqkit bedtools gawk; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool introuvable" >&2; exit 1; }
done
if command -v pigz >/dev/null; then GZ="pigz -p 8"; else GZ="gzip -3"; fi

RC=$(echo "$MOTIF" | tr 'ACGTRYKMBDHVNacgtrykmbdhvn' 'TGCAYRMKVHDBNtgcayrmkvhdbn' | rev)
if [[ "${RC^^}" != "${MOTIF^^}" ]]; then
  echo "ERROR: motif $MOTIF non palindromique (revcomp $RC) : non supporté" >&2; exit 1
fi
MOTIF_LEN=${#MOTIF}

POREC_FQ=$(realpath "$POREC_FQ"); HIFI_FQ=$(realpath "$HIFI_FQ")
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "======================================================================"
echo "      STARTING HIFI-GUIDED PORE-C DIGESTION PIPELINE (UNIX/BASH)"
echo "======================================================================"
echo "Pore-C File   : $POREC_FQ"
echo "HiFi File     : $HIFI_FQ"
echo "Motif         : $MOTIF (Cut Offset: $CUT_OFFSET, motif on both monomers: $DUP_MOTIF)"
echo "Protection    : >= $MIN_COV HiFi reads continuous over site +/- ${FLANK}bp (gaps < ${MAX_GAP}bp), otherwise cut"
echo "Min monomer   : ${MIN_MONO_LEN}bp (shorter monomers removed before pairing)"
echo "Threads       : $THREADS"
echo "Output dir    : $PWD"
echo "======================================================================"

# ------------------------------------------------------------------------------
# STEP 1: LONGUEURS DES READS PORE-C (ordre du FASTQ = ordre de tri pour bedtools)
# ------------------------------------------------------------------------------
if [[ ! -e step1.done ]]; then
  echo -e "\n[STEP 1/5] Extracting Pore-C read lengths..."
  seqkit fx2tab -n -i -l -j "$THREADS" "$POREC_FQ" | cut -f1,2 > porec.genome
  touch step1.done
fi
read -r N_POREC TOT_BP < <(awk '{n++; s += $2} END {print n + 0, s + 0}' porec.genome)
echo "  -> Pore-C reads input        : $N_POREC"
echo "  -> Total bases / mean length : $TOT_BP bp / $(awk -v t="$TOT_BP" -v n="$N_POREC" 'BEGIN {printf "%.1f", n ? t / n : 0}') bp"

# ------------------------------------------------------------------------------
# STEP 2: ALIGNEMENT PORE-C (requête) SUR LES READS HIFI (cible) -> BLOCS CONTINUS
# ------------------------------------------------------------------------------
# Chaque monomère Pore-C s'aligne sur les reads HiFi de son locus ; -N garde les
# autres reads HiFi du locus (secondaires). -c donne le CIGAR : un alignement est
# coupé en blocs à chaque indel >= MAX_GAP. Puis, pour un même read HiFi (même
# brin), les blocs colinéaires séparés de < MAX_GAP (sur le Pore-C ET sur le HiFi)
# sont fusionnés : un read HiFi aligné sur les deux flancs d'un site, avec un
# petit trou dû aux erreurs de séquençage, reste continu à travers ce site.
# Sortie (ordre du FASTQ, triée par début dans chaque read) :
#   read_porec  début  fin  read_hifi
BLOCKS=hifi_blocks.bed.gz
if [[ ! -e step2.done ]]; then
  echo -e "\n[STEP 2/5] Aligning Pore-C reads on HiFi reads (minimap2 -c, $MM2_PRESET)..."
  mkdir -p mm2_tmp
  minimap2 -c -x "$MM2_PRESET" -t "$THREADS" -I "$MM2_BATCH" --split-prefix mm2_tmp/split \
      --secondary=yes -N "$MM2_N" $MM2_EXTRA "$HIFI_FQ" "$POREC_FQ" 2> minimap2.log | \
  gawk -v G="$MAX_GAP" -v LONG="[0-9]{${#MAX_GAP},}[IDN]" '
    BEGIN { OFS = "\t" }
    function add(q1, q2, t1, t2) {
      if (q2 <= q1) return
      n++; HT[n] = hifi; HS[n] = st; Q1[n] = q1; Q2[n] = q2; T1[n] = t1; T2[n] = t2
    }
    function flush(   i, p, m, gq, gt) {
      if (!n) return
      # fusion des blocs colinéaires d un même read HiFi
      delete K
      for (i = 1; i <= n; i++) K[i] = sprintf("%s\t%s\t%012d", HT[i], HS[i], Q1[i])
      PROCINFO["sorted_in"] = "@val_str_asc"
      p = 0; m = 0
      for (i in K) {
        if (p && HT[i] == HT[p] && HS[i] == HS[p]) {
          gq = Q1[i] - Q2[p]
          gt = (HS[i] == "+") ? T1[i] - T2[p] : T1[p] - T2[i]
          if (gq > -G && gq < G && gt > -G && gt < G) {
            if (Q2[i] > Q2[p]) Q2[p] = Q2[i]
            if (HS[i] == "+") T2[p] = T2[i]; else T1[p] = T1[i]
            continue
          }
        }
        if (p) { m++; MQ1[m] = Q1[p]; MQ2[m] = Q2[p]; MT[m] = HT[p] }
        p = i
      }
      m++; MQ1[m] = Q1[p]; MQ2[m] = Q2[p]; MT[m] = HT[p]
      # sortie triée par début sur le read Pore-C (requis par bedtools map)
      delete O
      for (i = 1; i <= m; i++) O[i] = MQ1[i]
      PROCINFO["sorted_in"] = "@val_num_asc"
      for (i in O) print cur, MQ1[i], MQ2[i], MT[i]
      n = 0
    }
    $1 != cur { flush(); cur = $1 }
    {
      st = $5; hifi = $6
      cg = ""
      for (i = 13; i <= NF; i++) if (substr($i, 1, 5) == "cg:Z:") { cg = substr($i, 6); break }
      if (cg !~ LONG) { add($3 + 0, $4 + 0, $8 + 0, $9 + 0); next }   # aucun indel >= MAX_GAP
      k = split(cg, L, /[MIDNSHP=X]/, OP)
      qp = (st == "+") ? $3 : $4; tp = $8 + 0; bq = qp; bt = tp
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
    END { flush() }' | $GZ > "$BLOCKS"
  rm -rf mm2_tmp
  touch step2.done
fi
N_ALN_READS=$($GZ -dc "$BLOCKS" | cut -f1 | uniq | wc -l)
echo "  -> Pore-C reads with HiFi alignments : $N_ALN_READS"

# ------------------------------------------------------------------------------
# STEP 3: SITES DE RESTRICTION (position de coupure exacte)
# ------------------------------------------------------------------------------
# sites.tsv : read  longueur  début_motif  fin_motif  x   (x = coupure = début + CUT_OFFSET)
if [[ ! -e step3.done ]]; then
  echo -e "\n[STEP 3/5] Locating restriction motifs ($MOTIF)..."
  seqkit locate -P -d -i -j "$THREADS" -p "$MOTIF" --bed "$POREC_FQ" | \
  awk -v off="$CUT_OFFSET" -v G=porec.genome 'BEGIN { OFS = "\t" }
    {
      while ($1 != r) {
        if ((getline line < G) <= 0) { print "ERROR: read " $1 " absent/out of order in porec.genome" > "/dev/stderr"; exit 1 }
        split(line, a, "\t"); r = a[1]; len = a[2]
      }
      x = $2 + off
      if (x > 0 && x < len) print $1, len, $2, $3, x      # x = 0 ou x = longueur : rien à couper
    }' > sites.tsv
  touch step3.done
fi
N_MOTIFS=$(if [[ -e step4.done ]]; then $GZ -dc sites.final.tsv.gz | wc -l; else wc -l < sites.tsv; fi)
echo "  -> Candidate motif sites     : $N_MOTIFS"

# ------------------------------------------------------------------------------
# STEP 4: PROTECTION DES SITES PAR LES READS HIFI (bedtools map)
# ------------------------------------------------------------------------------
#  span = nb de reads HiFi DISTINCTS dont un bloc continu couvre ENTIÈREMENT
#         [x-FLANK, x+FLANK] (fenêtre tronquée aux extrémités du read)
#  P (protégé) : span >= MIN_COV -> pas de coupe
#  C (coupé)   : span <  MIN_COV -> coupe, y compris sans aucun read HiFi
# sites.final.tsv.gz : read  longueur  début_motif  fin_motif  x  span  classe
if [[ ! -e step4.done ]]; then
  echo -e "\n[STEP 4/5] Computing HiFi protection at each site (bedtools map)..."
  paste sites.tsv \
    <(awk -v f="$FLANK" 'BEGIN { OFS = "\t" } { s = $5 - f; e = $5 + f; print $1, (s < 0 ? 0 : s), (e > $2 ? $2 : e) }' sites.tsv \
        | bedtools map -g porec.genome -b "$BLOCKS" -c 4 -o count_distinct -f 1.0 -a - | cut -f4) | \
  awk -v mc="$MIN_COV" 'BEGIN { OFS = "\t" }
    NF != 6 { print "ERROR: bedtools map output truncated at line " NR > "/dev/stderr"; exit 1 }
    { print $0, ($6 >= mc ? "P" : "C") }' | $GZ > sites.final.tsv.gz
  [[ $($GZ -dc sites.final.tsv.gz | wc -l) -eq $N_MOTIFS ]] || { echo "ERROR: sites.final.tsv.gz incomplete" >&2; exit 1; }
  rm -f sites.tsv
  touch step4.done
fi

read -r N_PROT N_CUT < <($GZ -dc sites.final.tsv.gz | awk '$7 == "P" { p++ } $7 == "C" { c++ } END { print p + 0, c + 0 }')
pct() { awk -v a="$1" -v b="$2" 'BEGIN { if (b > 0) printf "%.2f", 100 * a / b; else print 0 }'; }
echo "  -> Protected sites (>= $MIN_COV HiFi) : $N_PROT ($(pct "$N_PROT" "$N_MOTIFS")%)"
echo "  -> Cut sites (< $MIN_COV HiFi)        : $N_CUT ($(pct "$N_CUT" "$N_MOTIFS")%)"

# ------------------------------------------------------------------------------
# STEP 5: DÉCOUPE + FILTRE DES PSEUDO-MONOMÈRES + PAIRES ALL-TO-ALL
# ------------------------------------------------------------------------------
# Tous les sites non protégés sont coupés. Un motif voisin d'une vraie jonction
# (à moins de ~FLANK pb) n'est pas protégé non plus : entre les deux coupes naît
# un pseudo-monomère de quelques pb. Les monomères < MIN_MONO_LEN sont écartés
# AVANT l'appariement (ils ne sont de toute façon pas mappables).
# Un read de n monomères retenus donne n(n-1)/2 paires @read:i-j/1 et /2.
if [[ ! -e step5.done ]]; then
  echo -e "\n[STEP 5/5] Cutting reads & writing all-to-all pseudo-Hi-C pairs..."
  seqkit fx2tab -i -j "$THREADS" "$POREC_FQ" | \
  awk -v CUTS=<($GZ -dc sites.final.tsv.gz | awk '$7 == "C" { print $1 "\t" $3 "\t" $4 "\t" $5 }') \
      -v dup="$DUP_MOTIF" -v minlen="$MIN_MONO_LEN" -v wm="$WRITE_MONOMERS" \
      -v c1="$GZ > ${PREFIX}_R1.fastq.gz" -v c2="$GZ > ${PREFIX}_R2.fastq.gz" -v cm="$GZ > ${PREFIX}_monomers.fastq.gz" '
    BEGIN { FS = "\t"; more = ((getline cl < CUTS) > 0) }
    {
      id = $1; seq = $2; qual = $3; L = length(seq); nc = 0
      while (more) {
        split(cl, a, "\t")
        if (a[1] != id) break
        nc++; X[nc] = a[4]; MS[nc] = a[2]; ME[nc] = a[3]
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
      print reads + 0, mono + 0, short + 0, kept + 0, multi + 0, single + 0, pairs + 0 > "step5.stats"
    }'
  touch step5.done
fi

read -r READS MONO SHORT KEPT MULTI SINGLE PAIRS < step5.stats
echo "  -> Monomers generated            : $MONO"
echo "  -> Monomers < ${MIN_MONO_LEN}bp removed     : $SHORT"
echo "  -> Monomers kept                 : $KEPT ($(awk -v m="$KEPT" -v r="$READS" 'BEGIN {printf "%.2f", r ? m / r : 0}') monomers/read)"
echo "  -> Multi-monomer reads           : $MULTI"
echo "  -> Reads with < 2 monomers       : $SINGLE (no pair possible)"
echo "  -> Pseudo-Hi-C pairs (all-to-all): $PAIRS"
echo "  -> Output files                  : $PWD/${PREFIX}_R1.fastq.gz ($(du -h "${PREFIX}_R1.fastq.gz" | cut -f1))"
echo "                                     $PWD/${PREFIX}_R2.fastq.gz ($(du -h "${PREFIX}_R2.fastq.gz" | cut -f1))"
echo "  -> Per-site table                : $PWD/sites.final.tsv.gz"
echo "     (read len motif_start motif_end cut span class)"
echo "======================================================================"
echo "                   ALL PIPELINE STEPS COMPLETED"
echo "======================================================================"
