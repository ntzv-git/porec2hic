#!/usr/bin/env bash
# ==============================================================================
#  HIFI-GUIDED PORE-C DIGESTION -> ALL-TO-ALL PSEUDO-HI-C PAIRS
# ==============================================================================
#  Usage : POREC_FILE=porec.fq.gz HIFI_FQ=hifi.fq.gz ./porec2hic_hifi.sh
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
POREC_FILE=${POREC_FILE:?"POREC_FILE non défini"}
HIFI_FQ=${HIFI_FQ:?"HIFI_FQ non défini (FASTQ/FASTA HiFi ou index .mmi)"}
OUTDIR=${OUTDIR:-porec2hic_out}
PREFIX=${PREFIX:-porec_hic}

MOTIF=${MOTIF:-CATG}              # motif palindromique, IUPAC accepté (GATC, GANTC, ...)
CUT_OFFSET=${CUT_OFFSET:-4}       # coupure dans le motif : NlaIII CATG^ = 4, DpnII ^GATC = 0
DUP_MOTIF=${DUP_MOTIF:-1}         # 1 = le motif reconstitué à la ligation est gardé sur les 2 monomères

MIN_COV=${MIN_COV:-3}             # reads HiFi continus sur [x-FLANK, x+FLANK] pour PROTÉGER le site
FLANK=${FLANK:-15}                # marge (pb) exigée de part et d'autre du site (> longueur motif + débordement)
MAX_GAP=${MAX_GAP:-50}            # indel/trou (pb) au-delà duquel un alignement HiFi n'est plus continu
MIN_SIDE_COV=${MIN_SIDE_COV:-$MIN_COV}  # reads HiFi d'un côté du site pour valider une JONCTION
SIDE_WINDOW=${SIDE_WINDOW:-100}   # fenêtre (pb) de chaque côté du site pour compter ces reads
CUT_UNRESOLVED=${CUT_UNRESOLVED:-0}     # 1 = couper aussi les sites sans information HiFi
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

POREC_FILE=$(realpath "$POREC_FILE"); HIFI_FQ=$(realpath "$HIFI_FQ")
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "======================================================================"
echo "      STARTING HIFI-GUIDED PORE-C DIGESTION PIPELINE (UNIX/BASH)"
echo "======================================================================"
echo "Pore-C File   : $POREC_FILE"
echo "HiFi File     : $HIFI_FQ"
echo "Motif         : $MOTIF (Cut Offset: $CUT_OFFSET, motif on both monomers: $DUP_MOTIF)"
echo "Protection    : >= $MIN_COV HiFi reads continuous over site +/- ${FLANK}bp (gaps < ${MAX_GAP}bp)"
echo "Junction      : < $MIN_COV continuous and >= $MIN_SIDE_COV HiFi reads within ${SIDE_WINDOW}bp on one side"
echo "Threads       : $THREADS"
echo "Output dir    : $PWD"
echo "======================================================================"

# ------------------------------------------------------------------------------
# STEP 1: LONGUEURS DES READS PORE-C (ordre du FASTQ = ordre de tri pour bedtools)
# ------------------------------------------------------------------------------
if [[ ! -e step1.done ]]; then
  echo -e "\n[STEP 1/6] Extracting Pore-C read lengths..."
  seqkit fx2tab -n -i -l -j "$THREADS" "$POREC_FILE" | cut -f1,2 > porec.genome
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
  echo -e "\n[STEP 2/6] Aligning Pore-C reads on HiFi reads (minimap2 -c, $MM2_PRESET)..."
  mkdir -p mm2_tmp
  minimap2 -c -x "$MM2_PRESET" -t "$THREADS" -I "$MM2_BATCH" --split-prefix mm2_tmp/split \
      --secondary=yes -N "$MM2_N" $MM2_EXTRA "$HIFI_FQ" "$POREC_FILE" 2> minimap2.log | \
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
  echo -e "\n[STEP 3/6] Locating restriction motifs ($MOTIF)..."
  seqkit locate -P -d -i -j "$THREADS" -p "$MOTIF" --bed "$POREC_FILE" | \
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
N_MOTIFS=$(if [[ -e step5.done ]]; then $GZ -dc sites.final.tsv.gz | wc -l; else wc -l < sites.tsv; fi)
echo "  -> Candidate motif sites     : $N_MOTIFS"

# ------------------------------------------------------------------------------
# STEP 4: CLASSIFICATION DES SITES PAR LES READS HIFI (bedtools map)
# ------------------------------------------------------------------------------
#  span = reads HiFi distincts dont un bloc continu couvre ENTIÈREMENT [x-FLANK, x+FLANK]
#  covL = reads HiFi distincts avec un bloc dans [x-SIDE_WINDOW, x)
#  covR = reads HiFi distincts avec un bloc dans [x, x+SIDE_WINDOW)
#  P (protégé)   : span >= MIN_COV
#  J (jonction)  : span <  MIN_COV et max(covL, covR) >= MIN_SIDE_COV
#  U (non résolu): sinon (pas assez de HiFi autour du site)
MAP="bedtools map -g porec.genome -b $BLOCKS -c 4 -o count_distinct"
if [[ ! -e step4.done ]]; then
  echo -e "\n[STEP 4/6] Computing HiFi continuity at each site (bedtools map)..."
  paste sites.tsv \
    <(awk -v f="$FLANK" 'BEGIN { OFS = "\t" } { s = $5 - f; e = $5 + f; print $1, (s < 0 ? 0 : s), (e > $2 ? $2 : e) }' sites.tsv \
        | $MAP -f 1.0 -a - | cut -f4) \
    <(awk -v w="$SIDE_WINDOW" 'BEGIN { OFS = "\t" } { s = $5 - w; print $1, (s < 0 ? 0 : s), $5 }' sites.tsv \
        | $MAP -a - | cut -f4) \
    <(awk -v w="$SIDE_WINDOW" 'BEGIN { OFS = "\t" } { e = $5 + w; print $1, $5, (e > $2 ? $2 : e) }' sites.tsv \
        | $MAP -a - | cut -f4) | \
  awk -v mc="$MIN_COV" -v ms="$MIN_SIDE_COV" 'BEGIN { OFS = "\t" }
    NF != 8 { print "ERROR: bedtools map output truncated at line " NR > "/dev/stderr"; exit 1 }
    { c = ($6 >= mc) ? "P" : (($7 >= ms || $8 >= ms) ? "J" : "U"); print $0, c }' > sites.cls.tsv
  [[ $(wc -l < sites.cls.tsv) -eq $N_MOTIFS ]] || { echo "ERROR: sites.cls.tsv incomplete" >&2; exit 1; }
  touch step4.done
fi

# ------------------------------------------------------------------------------
# STEP 5: UN SEUL SITE PAR JONCTION
# ------------------------------------------------------------------------------
# Un autre motif à moins de 2*FLANK d'une jonction n'est pas non plus traversé
# par les HiFi. Dans une grappe de sites J rapprochés, on garde le site le plus
# proche du point de cassure = milieu entre la médiane des fins des blocs HiFi
# à gauche et la médiane des débuts des blocs HiFi à droite. Les autres -> "j".
if [[ ! -e step5.done ]]; then
  echo -e "\n[STEP 5/6] Resolving clusters of junction sites..."
  awk -v f="$FLANK" 'BEGIN { OFS = "\t" }
    function out(   i) { if (nb > 1) for (i = 1; i <= nb; i++) print buf[i]; nb = 0 }
    $9 != "J" { out(); next }
    { if (nb && ($1 != pr || $5 - px >= 2 * f)) out(); buf[++nb] = $1 "\t" $2 "\t" $5; pr = $1; px = $5 }
    END { out() }' sites.cls.tsv > clusters.tsv
  # fins de blocs autour de x (gauche) et débuts de blocs autour de x (droite)
  paste clusters.tsv \
    <(awk -v f="$FLANK" -v w="$SIDE_WINDOW" 'BEGIN { OFS = "\t" } { s = $3 - w; e = $3 + f; print $1, (s < 0 ? 0 : s), (e > $2 ? $2 : e) }' clusters.tsv \
        | bedtools map -g porec.genome -b "$BLOCKS" -c 3 -o collapse -a - | cut -f4) \
    <(awk -v f="$FLANK" -v w="$SIDE_WINDOW" 'BEGIN { OFS = "\t" } { s = $3 - f; e = $3 + w; print $1, (s < 0 ? 0 : s), (e > $2 ? $2 : e) }' clusters.tsv \
        | bedtools map -g porec.genome -b "$BLOCKS" -c 2 -o collapse -a - | cut -f4) | \
  gawk -v f="$FLANK" -v w="$SIDE_WINDOW" 'BEGIN { OFS = "\t" }
    function med(list, lo, hi,   n, i, v, k) {
      n = split(list, v, ","); k = 0; delete M
      for (i = 1; i <= n; i++) if (v[i] != "." && v[i] >= lo && v[i] <= hi) M[++k] = v[i] + 0
      if (!k) return ""
      asort(M)
      return M[int((k + 1) / 2)]
    }
    {
      e = med($4, $3 - w, $3 + f); s = med($5, $3 - f, $3 + w)
      b = (e != "" && s != "") ? (e + s) / 2 : (e != "" ? e : (s != "" ? s : $3))
      print $1, $3, b
    }' > clusters.bp.tsv

  awk -v f="$FLANK" -v BP=clusters.bp.tsv 'BEGIN { OFS = "\t"; while ((getline l < BP) > 0) { split(l, a, "\t"); B[a[1] SUBSEP a[2]] = a[3] } }
    function d(i,   v) { v = xs[i] - B[rd SUBSEP xs[i]]; return v < 0 ? -v : v }
    function out(   i, best) {
      if (nb > 1) { best = 1; for (i = 2; i <= nb; i++) if (d(i) < d(best)) best = i
                    for (i = 1; i <= nb; i++) if (i != best) sub(/J$/, "j", buf[i]) }
      for (i = 1; i <= nb; i++) print buf[i]
      nb = 0
    }
    $9 != "J" { out(); print; next }
    { if (nb && ($1 != rd || $5 - xs[nb] >= 2 * f)) out(); rd = $1; buf[++nb] = $0; xs[nb] = $5 }
    END { out() }' sites.cls.tsv | $GZ > sites.final.tsv.gz
  rm -f clusters.tsv clusters.bp.tsv sites.cls.tsv sites.tsv
  touch step5.done
fi

declare -A C=([P]=0 [J]=0 [j]=0 [U]=0)
while read -r k v; do C[$k]=$v; done < <($GZ -dc sites.final.tsv.gz | awk '{c[$9]++} END {for (k in c) print k, c[k]}')
pct() { awk -v a="$1" -v b="$2" 'BEGIN { if (b > 0) printf "%.2f", 100 * a / b; else print 0 }'; }
echo "  -> Protected sites (P)          : ${C[P]} ($(pct "${C[P]}" "$N_MOTIFS")%)"
echo "  -> HiFi-validated junctions (J) : ${C[J]} ($(pct "${C[J]}" "$N_MOTIFS")%)"
echo "  -> Motifs next to a junction (j): ${C[j]} ($(pct "${C[j]}" "$N_MOTIFS")%)  -> not cut"
echo "  -> Unresolved, no HiFi info (U) : ${C[U]} ($(pct "${C[U]}" "$N_MOTIFS")%)  -> $([[ $CUT_UNRESOLVED == 1 ]] && echo cut || echo 'not cut')"

# ------------------------------------------------------------------------------
# STEP 6: DÉCOUPE + PAIRES ALL-TO-ALL (seqkit fx2tab + awk, une passe, sans filtre)
# ------------------------------------------------------------------------------
# Tous les monomères sont conservés, quelle que soit leur taille : le filtrage se
# fait en aval sur les paires. Un read de n monomères donne n(n-1)/2 paires.
if [[ ! -e step6.done ]]; then
  echo -e "\n[STEP 6/6] Cutting reads & writing all-to-all pseudo-Hi-C pairs..."
  CUT_CLASSES="J"; [[ "$CUT_UNRESOLVED" == 1 ]] && CUT_CLASSES="JU"
  seqkit fx2tab -i -j "$THREADS" "$POREC_FILE" | \
  awk -v CUTS=<($GZ -dc sites.final.tsv.gz | awk -v cc="$CUT_CLASSES" 'index(cc, $9) { print $1 "\t" $3 "\t" $4 "\t" $5 }') \
      -v dup="$DUP_MOTIF" -v wm="$WRITE_MONOMERS" \
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
      for (k = 1; k <= nc; k++) {
        end = dup ? ME[k] : X[k]
        if (end > beg) { n++; S[n] = beg; E[n] = end }
        beg = dup ? MS[k] : X[k]
      }
      if (L > beg) { n++; S[n] = beg; E[n] = L }
      reads++; mono += n
      if (n == 1) single++
      for (i = 1; i <= n; i++) {
        sq[i] = substr(seq, S[i] + 1, E[i] - S[i]); ql[i] = substr(qual, S[i] + 1, E[i] - S[i])
        if (wm) print "@" id "_" S[i] "_" E[i] "\n" sq[i] "\n+\n" ql[i] | cm
      }
      if (n >= 2) {
        multi++
        for (i = 1; i <= n; i++)
          for (j = i + 1; j <= n; j++) {
            tag = "@" id ":" i "-" j
            print tag "/1\n" sq[i] "\n+\n" ql[i] | c1
            print tag "/2\n" sq[j] "\n+\n" ql[j] | c2
            pairs++
          }
      }
    }
    END {
      if (more) { print "ERROR: cut list not exhausted (" cl ") -> order mismatch with the FASTQ" > "/dev/stderr"; exit 1 }
      close(c1); close(c2); if (wm) close(cm)
      print reads + 0, mono + 0, multi + 0, single + 0, pairs + 0 > "step6.stats"
    }'
  touch step6.done
fi

read -r READS MONO MULTI SINGLE PAIRS < step6.stats
echo "  -> Monomers generated          : $MONO ($(awk -v m="$MONO" -v r="$READS" 'BEGIN {printf "%.2f", r ? m / r : 0}') monomers/read)"
echo "  -> Multi-monomer reads         : $MULTI"
echo "  -> Single-monomer reads        : $SINGLE (no pair possible)"
echo "  -> Pseudo-Hi-C pairs (all-to-all): $PAIRS"
echo "  -> Output files                : $PWD/${PREFIX}_R1.fastq.gz ($(du -h "${PREFIX}_R1.fastq.gz" | cut -f1))"
echo "                                   $PWD/${PREFIX}_R2.fastq.gz ($(du -h "${PREFIX}_R2.fastq.gz" | cut -f1))"
echo "  -> Per-site table              : $PWD/sites.final.tsv.gz"
echo "     (read len motif_start motif_end cut span covL covR class)"
echo "======================================================================"
echo "                   ALL PIPELINE STEPS COMPLETED"
echo "======================================================================"
