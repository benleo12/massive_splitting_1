#!/bin/bash
# Collinear-limit scan: fix s,u and m^2, drive t -> 0, and evaluate the masters at
# each point (fast path, runtime kinematics).  Used to show that a nonzero quark mass
# regulates the t->0 (collinear) region, whereas the massless amplitude diverges.
#   usage: ./collinear_scan.sh <mq2> <outdir>
set -e
MQ2=${1:-1}; OUT=${2:-/tmp/scan_mq$1}
cd "$(dirname "$0")"
mkdir -p "$OUT"
S=-3; U=-2
IDS="B0_t_0_mq B0_u_0_mq B0_mq2_0_mq B0_s_mq_mq B0_Q_mq_mq C0_0_mq_t C0_0_mq_u C0_0_s_Q C0_mq_mq_s C0_mq_t_Q C0_mq_u_Q D0_a D0_b"
for T in -1 -0.3 -0.1 -0.03 -0.01 -0.003; do
  f="$OUT/t${T}.txt"; : > "$f"
  for id in $IDS; do
    docker run --rm -v "$PWD":/work -w /work pysecdec:1mass \
      python3 psd_fast.py eval $id $S $T $U $MQ2 2>/dev/null | grep PSDRESULT >> "$f"
  done
  n=$(grep -c PSDRESULT "$f")
  python3 parse_psd.py "$OUT/m_t${T}.m" "$f" >/dev/null
  echo "  t=$T  ($n/13 masters)"
done
echo "SCAN_DONE mq2=$MQ2 -> $OUT"
