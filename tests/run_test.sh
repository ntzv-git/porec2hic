#!/usr/bin/env bash
# End-to-end test on simulated data (needs: pip install mappy).
# minimap2 is replaced by a mappy-based stand-in; FAKE_CHAIN_JITTER=30 mimics
# the imprecise ends of chain-only PAF (minimap2 without -c).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
W=${1:-$(mktemp -d)}
mkdir -p "$W/bin"
install -m 755 "$HERE/fake_minimap2.py" "$W/bin/minimap2"
cd "$W"
python3 "$HERE/simulate.py" "$W"
FAKE_CHAIN_JITTER=${FAKE_CHAIN_JITTER:-30} PATH="$W/bin:$PATH" THREADS=4 MM2_N=150 WRITE_SITES=1 \
  POREC_FILE=porec.fq.gz HIFI_FQ=hifi.fq.gz OUTDIR=out bash "$HERE/../porec2hic_hifi.sh"
python3 "$HERE/eval.py"
