#!/bin/bash
# FAST evaluation of all 15 masters at one point, using pre-built libraries
# (run ./build_masters_2mass.sh ONCE first).  ~sub-second/point native, no recompile.
# Uses native pySecDec if importable, else the Docker image.
#   usage:  ./run_masters_fast_2mass.sh <s> <t> <u> <m1sq> <m2sq> [outfile.m]
set -e
S=$1; T=$2; U=$3; M1=$4; M2=$5; OUT=${6:-masters_2mass.m}
cd "$(dirname "$0")"
if [ ! -e m_D0_b/m_D0_b_pylink.so ]; then
  echo "libraries not built yet -- run ./build_masters_2mass.sh first"; exit 1
fi
if python3 -c "import pySecDec" 2>/dev/null; then NATIVE=1; fi
psd() { if [ -n "$NATIVE" ]; then python3 "$@";
        else docker run --rm -v "$PWD":/work -w /work pysecdec:1mass python3 "$@"; fi; }
RAW=$(mktemp)
psd psd_fast_2mass.py evalall $S $T $U $M1 $M2 2>/dev/null | grep PSDRESULT > "$RAW"
python3 parse_psd.py "$OUT" "$RAW"
rm -f "$RAW"
echo "[fast] wrote $OUT"
