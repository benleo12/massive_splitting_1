#!/bin/bash
# run all 13 masters of the 1-mass squark amplitude via pySecDec, at one point,
# recompiling per point (simple path).  Uses native pySecDec if importable, else Docker.
#   usage:  ./run_masters.sh <s> <t> <u> <mq2> [outfile.m]
# Euclidean points (s,t,u<0, mq2>0) are safest.
set -e
S=$1; T=$2; U=$3; MQ2=$4; OUT=${5:-masters.m}
cd "$(dirname "$0")"
if python3 -c "import pySecDec" 2>/dev/null; then NATIVE=1; fi
psd() { if [ -n "$NATIVE" ]; then python3 "$@";
        else docker run --rm -v "$PWD":/work -w /work pysecdec:1mass python3 "$@"; fi; }
RAW=$(mktemp)
IDS="B0_t_0_mq B0_u_0_mq B0_mq2_0_mq B0_s_mq_mq B0_Q_mq_mq \
     C0_0_mq_t C0_0_mq_u C0_0_s_Q C0_mq_mq_s C0_mq_t_Q C0_mq_u_Q D0_a D0_b"
for id in $IDS; do
  echo "[run] $id  (s=$S t=$T u=$U mq2=$MQ2)"
  psd psd_master.py $id $S $T $U $MQ2 2>/dev/null | grep PSDRESULT >> "$RAW"
done
python3 parse_psd.py "$OUT" "$RAW"
rm -f "$RAW"
echo "[run] wrote $OUT"
