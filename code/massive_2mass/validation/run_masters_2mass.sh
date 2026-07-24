#!/bin/bash
# run all 15 masters of the 2-mass squark amplitude via pySecDec, at one point,
# recompiling per point (simple path).  Uses native pySecDec if importable, else Docker.
#   usage:  ./run_masters_2mass.sh <s> <t> <u> <m1sq> <m2sq> [outfile.m]
# Euclidean points (s,t,u<0, m1sq,m2sq>0) are safest.
set -e
S=$1; T=$2; U=$3; M1=$4; M2=$5; OUT=${6:-masters_2mass.m}
cd "$(dirname "$0")"
if python3 -c "import pySecDec" 2>/dev/null; then NATIVE=1; fi
psd() { if [ -n "$NATIVE" ]; then python3 "$@";
        else docker run --rm -v "$PWD":/work -w /work pysecdec:1mass python3 "$@"; fi; }
RAW=$(mktemp)
IDS="B0_m1_0_m1 B0_m2_0_m2 B0_s_m1_m2 B0_Q_m1_m2 B0_t_0_m1 B0_u_0_m2 \
     C0_0_m1_t C0_0_m2_u C0_0_s_Q_a C0_0_s_Q_b C0_m1_m2_s C0_m1_u_Q C0_m2_t_Q \
     D0_a D0_b"
for id in $IDS; do
  echo "[run] $id  (s=$S t=$T u=$U m1sq=$M1 m2sq=$M2)"
  psd psd_master_2mass.py $id $S $T $U $M1 $M2 2>/dev/null | grep PSDRESULT >> "$RAW"
done
python3 parse_psd.py "$OUT" "$RAW"
rm -f "$RAW"
echo "[run] wrote $OUT"
