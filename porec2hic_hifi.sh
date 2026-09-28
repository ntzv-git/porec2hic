#!/usr/bin/env bash
# ==============================================================================
#  HIFI-GUIDED PORE-C DIGESTION -> ALL-TO-ALL PSEUDO-HI-C PAIRS
# ==============================================================================
#  Usage : POREC_FILE=porec.fq.gz HIFI_FQ=hifi.fq.gz ./porec2hic_hifi.sh
#  Toutes les variables ci-dessous peuvent être surchargées par l'environnement.
#  Chaque étape est reprise si son fichier .done existe (reprise après crash).
# ==============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# VARIABLES DU PIPELINE
# ------------------------------------------------------------------------------
POREC_FILE=${POREC_FILE:?"POREC_FILE non défini"}
HIFI_FQ=${HIFI_FQ:?"HIFI_FQ non défini (FASTQ/FASTA HiFi ou index .mmi)"}
OUTDIR=${OUTDIR:-porec2hic_out}
PREFIX=${PREFIX:-porec_hic}

MOTIF=${MOTIF:-CATG}              # IUPAC accepté (ex: GATC, GANTC)
CUT_OFFSET=${CUT_OFFSET:-4}       # position de coupure dans le motif : NlaIII CATG^ = 4, DpnII ^GATC = 0
MIN_COV=${MIN_COV:-3}             # nb d'alignements HiFi traversant le site pour le protéger
MIN_SIDE_COV=${MIN_SIDE_COV:-$MIN_COV}  # nb d'alignements HiFi d'un côté pour valider une jonction
FLANK=${FLANK:-25}                # un alignement doit dépasser le site de FLANK pb de chaque côté
SIDE_WINDOW=${SIDE_WINDOW:-100}   # fenêtre (pb) de part et d'autre du site où compter les alignements HiFi
CUT_UNRESOLVED=${CUT_UNRESOLVED:-0}     # 1 = couper aussi les sites sans info HiFi (ancien comportement)
DUP_MOTIF=${DUP_MOTIF:-1}         # 1 = le motif (partagé à la ligation) est gardé sur les deux monomères
MIN_MONO_LEN=${MIN_MONO_LEN:-50}  # monomères plus courts ignorés (non mappables)
MAX_MONOMERS=${MAX_MONOMERS:-0}   # 0 = pas de limite ; sinon ignore les reads plus fragmentés
WRITE_MONOMERS=${WRITE_MONOMERS:-0}
WRITE_SITES=${WRITE_SITES:-0}     # table par site (span/covL/covR/classe) pour calibrer MIN_COV

THREADS=${THREADS:-96}
MM2_PRESET=${MM2_PRESET:-map-ont} # lr:hq pour ONT R10 Q20+ (minimap2 >= 2.27)
MM2_N=${MM2_N:-100}               # secondaires (= autres reads HiFi du même locus). Plafond GLOBAL par
                                  # read Pore-C (pas par monomère) : viser ~ MIN_COV x nb de monomères x 5
MM2_BATCH=${MM2_BATCH:-50G}       # taille des lots d'index HiFi (-I) : mémoire ~ 2-2.5 o/base
MM2_EXTRA=${MM2_EXTRA:-}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
POREC_FILE=$(realpath "$POREC_FILE"); HIFI_FQ=$(realpath "$HIFI_FQ")
mkdir -p "$OUTDIR"
cd "$OUTDIR"
AWK=$(command -v mawk || command -v awk)
if command -v pigz >/dev/null; then GZ="pigz -p 8"; else GZ="gzip -3"; fi

echo "======================================================================"
echo "      STARTING HIFI-GUIDED PORE-C DIGESTION PIPELINE (UNIX/BASH)"
echo "======================================================================"
echo "Pore-C File   : $POREC_FILE"
echo "HiFi File     : $HIFI_FQ"
echo "Motif         : $MOTIF (Cut Offset: $CUT_OFFSET, dup motif: $DUP_MOTIF)"
echo "Min Coverage  : span >= $MIN_COV -> protected ; side >= $MIN_SIDE_COV -> junction (flank ${FLANK}bp)"
echo "Threads       : $THREADS"
echo "Output dir    : $PWD"
echo "======================================================================"

# ------------------------------------------------------------------------------
# STEP 1: ALIGNEMENT PORE-C (requête) SUR LES READS HIFI (cible)
# ------------------------------------------------------------------------------
# Sens inversé par rapport à la version précédente : chaque monomère Pore-C est
# aligné sur les reads HiFi du même locus. minimap2 découpe naturellement
# l'alignement aux jonctions de ligation (primaire + supplémentaires), et -N
# garde les autres reads HiFi du locus (secondaires) => couverture HiFi par site
# sans dépendre du fait qu'un read HiFi "choisisse" ce read Pore-C.
# Sortie : PAF (sans CIGAR, pas de BAM/tri/index) compacté à 1 ligne par read :
#   read_id  qlen  qstart1,qstart2,...  qend1,qend2,...
# --split-prefix : index HiFi en lots de MM2_BATCH avec fusion correcte des hits.
ALN=hifi_aln.tsv.gz
if [[ ! -e step1.done ]]; then
  echo -e "\n[STEP 1/3] Aligning Pore-C reads on HiFi reads with minimap2 ($MM2_PRESET)..."
  mkdir -p mm2_tmp
  minimap2 -x "$MM2_PRESET" -t "$THREADS" -I "$MM2_BATCH" --split-prefix mm2_tmp/split \
      --secondary=yes -N "$MM2_N" $MM2_EXTRA "$HIFI_FQ" "$POREC_FILE" 2> minimap2.log | \
    $AWK 'BEGIN { OFS = "\t" }
      $1 != q { if (q != "") print q, ql, s, e; q = $1; ql = $2; s = $3; e = $4; next }
              { s = s "," $3; e = e "," $4 }
      END     { if (q != "") print q, ql, s, e }' | \
    $GZ > "$ALN"
  rm -rf mm2_tmp
  touch step1.done
fi
N_ALN_READS=$($GZ -dc "$ALN" | wc -l)
echo "  -> Pore-C reads with HiFi alignments : $N_ALN_READS"

# ------------------------------------------------------------------------------
# STEP 2: SITES + CLASSIFICATION HIFI + DÉCOUPE + PAIRES ALL-TO-ALL (1 passe)
# ------------------------------------------------------------------------------
if [[ ! -e step2.done ]]; then
  echo -e "\n[STEP 2/3] Classifying $MOTIF sites, cutting junctions & building all-to-all pairs..."
  OPTS=()
  [[ "$CUT_UNRESOLVED" == 1 ]] && OPTS+=(--cut-unresolved)
  [[ "$DUP_MOTIF" == 0 ]]      && OPTS+=(--no-dup-motif)
  [[ "$WRITE_MONOMERS" == 1 ]] && OPTS+=(--monomers)
  [[ "$WRITE_SITES" == 1 ]]    && OPTS+=(--sites)
  python3 "$SCRIPT_DIR/porec_hifi_split.py" \
      --fastq "$POREC_FILE" --aln "$ALN" --prefix "$PREFIX" \
      --motif "$MOTIF" --cut-offset "$CUT_OFFSET" \
      --min-cov "$MIN_COV" --min-side-cov "$MIN_SIDE_COV" --flank "$FLANK" --side-window "$SIDE_WINDOW" \
      --min-monomer-len "$MIN_MONO_LEN" --max-monomers "$MAX_MONOMERS" \
      -t "$THREADS" ${OPTS[@]+"${OPTS[@]}"}
  touch step2.done
fi

# ------------------------------------------------------------------------------
# STEP 3: RAPPORT
# ------------------------------------------------------------------------------
echo -e "\n[STEP 3/3] Summary"
declare -A S
while IFS=$'\t' read -r k v; do S[$k]=$v; done < "$PREFIX.stats.tsv"
pct() { $AWK -v a="$1" -v b="$2" 'BEGIN { if (b > 0) printf "%.2f", 100 * a / b; else print 0 }'; }
ratio() { $AWK -v a="$1" -v b="$2" 'BEGIN { if (b > 0) printf "%.2f", a / b; else print 0 }'; }
EVAL=$(( S[sites_total] - S[sites_edge] ))

echo "  -> Pore-C reads input            : ${S[reads]} ($(pct "${S[reads_with_aln]}" "${S[reads]}")% with HiFi alignments)"
echo "  -> Total bases / mean length     : ${S[bases]} bp / $(ratio "${S[bases]}" "${S[reads]}") bp"
echo "  -> Candidate $MOTIF sites         : ${S[sites_total]} (${S[sites_edge]} at read ends, ignored)"
echo "  -> Protected sites (span>=$MIN_COV)   : ${S[sites_protected]} ($(pct "${S[sites_protected]}" "$EVAL")%)"
echo "  -> HiFi-validated junctions      : ${S[sites_junction]} ($(pct "${S[sites_junction]}" "$EVAL")%)"
echo "  -> Junction-adjacent motifs (<2x$FLANK bp, not cut): ${S[sites_junction_adjacent]}"
echo "  -> Unresolved (no HiFi info)     : ${S[sites_unresolved]} ($(pct "${S[sites_unresolved]}" "$EVAL")%) $([[ $CUT_UNRESOLVED == 1 ]] && echo '-> cut' || echo '-> kept')"
echo "  -> Cuts applied                  : ${S[cuts]}"
echo "  -> Monomers (raw / <${MIN_MONO_LEN}bp / kept) : ${S[monomers_raw]} / ${S[monomers_short]} / ${S[monomers_kept]}"
echo "  -> Fragmentation ratio           : $(ratio "${S[monomers_kept]}" "${S[reads]}") monomers/read"
echo "  -> Multi-monomer reads           : ${S[reads_multi]}"
echo "  -> Single-monomer reads (dropped): ${S[reads_single]}  | empty: ${S[reads_empty]}  | capped: ${S[reads_capped]}"
echo "  -> Pseudo-Hi-C pairs (all-to-all): ${S[pairs]}"
echo "  -> Output files                  : $PWD/${PREFIX}_R1.fastq.gz ($(du -h "${PREFIX}_R1.fastq.gz" | cut -f1))"
echo "                                     $PWD/${PREFIX}_R2.fastq.gz ($(du -h "${PREFIX}_R2.fastq.gz" | cut -f1))"
echo "======================================================================"
echo "                   ALL PIPELINE STEPS COMPLETED"
echo "======================================================================"
