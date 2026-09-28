#!/usr/bin/env bash
# Test de bout en bout sur données simulées (hors pipeline : la simulation et
# l'évaluation utilisent python3, le pipeline lui-même non).
# Requiert minimap2, seqkit, bedtools, gawk dans le PATH.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
W=${1:-$(mktemp -d)}
mkdir -p "$W"
cd "$W"
[[ -s porec.fq.gz ]] || python3 "$HERE/simulate.py" "$W"
THREADS=${THREADS:-4} POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz OUTDIR=${OUTDIR:-out} \
  bash "$HERE/../porec2hic_hifi.sh"
python3 "$HERE/eval.py" "${OUTDIR:-out}"
