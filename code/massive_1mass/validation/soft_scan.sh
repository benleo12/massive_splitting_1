#!/bin/bash
# Soft-limit scan for the BCM (Bierenbaum-Czakon-Mitov, arXiv:1107.4384) cross-check.
#
# We fix the two massive emitters k1,k3 exactly on shell and scale the emitted gluon
# k2 = w*n  (n light-like), so that
#     s   = (k1+k3)^2 = 2 mq^2 + 2 pipj      is INDEPENDENT of w,
#     t   = mq^2 + 2 piq,   piq = w (k1.n)   -> 0 linearly,
#     u   = mq^2 + 2 pjq,   pjq = w (k3.n)   -> 0 linearly.
# This is PHYSICAL (Minkowski) kinematics -- pySecDec uses contour deformation and
# returns complex master values, which is what the comparison with BCM needs since
# their g_ij carries explicit absorptive (i*pi) parts.
#
#   usage: ./soft_scan.sh <outdir>
set -e
cd "$(dirname "$0")"
OUT=${1:-/tmp/softscan}
# how many CPUs the pySecDec containers may use (default 2; override: DCPUS=4 ./soft_scan.sh)
DCPUS=${DCPUS:-2}
mkdir -p "$OUT"

# the 14 masters of the non-abelian radiator (paveToIdN in amplitude_functions_nab.wl)
IDS="B0_mq2_0_mq B0_Q_mq_mq C0_ggt C0_0_mq_t C0_ggu C0_0_mq_u C0_0_s_Q \
     C0_mq_t_Q C0_mq_u_Q D0_a D0_b D0_c D0_d D0_e"

# point A: P=2, cos(theta)=1/2  -> s=20 ; point B: P=3, cos(theta)=-3/10 -> s=40
# fields: tag s mq2 dt du   with  t = mq2 + dt*w ,  u = mq2 + du*w
POINTS="C:100:1:7.5505102572168215:12.449489742783178 D:8:1:1.8284271247461903:3.8284271247461903"
OMEGAS="0.1 0.05 0.02 0.01"

for pt in $POINTS; do
  tag=$(echo $pt|cut -d: -f1); S=$(echo $pt|cut -d: -f2); MQ2=$(echo $pt|cut -d: -f3)
  DT=$(echo $pt|cut -d: -f4);  DU=$(echo $pt|cut -d: -f5)
  for W in $OMEGAS; do
    T=$(python3 -c "print(repr($MQ2 + $DT*$W))")
    U=$(python3 -c "print(repr($MQ2 + $DU*$W))")
    f="$OUT/${tag}_w${W}.txt"; : > "$f"
    for id in $IDS; do
      # CPU-capped: pySecDec/Qmc will otherwise saturate every core and make the
      # machine unusable.  --cpus caps the container at the cgroup level, which is
      # what actually protects interactive responsiveness; nice lowers priority too.
      docker run --rm --cpus="$DCPUS" --memory=4g -v "$PWD":/work -w /work pysecdec:1mass \
        nice -n 15 python3 psd_fast.py eval $id $S $T $U $MQ2 2>/dev/null | grep PSDRESULT >> "$f" || true
    done
    n=$(grep -c PSDRESULT "$f" || true)
    python3 parse_psd_cplx.py "$OUT/m_${tag}_w${W}.m" "$f" >/dev/null
    echo "  $tag  w=$W  s=$S t=$T u=$U   ($n/14 masters)"
  done
done
echo "SOFTSCAN_DONE -> $OUT"
