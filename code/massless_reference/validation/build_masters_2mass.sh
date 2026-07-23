#!/bin/bash
# ONE-TIME build of all 15 master libraries with kinematics as runtime parameters.
# Uses native pySecDec if importable (pip install pySecDec), else the Docker image.
#   usage:  ./build_masters_2mass.sh
set -e
cd "$(dirname "$0")"
if python3 -c "import pySecDec" 2>/dev/null; then NATIVE=1; fi
psd() { if [ -n "$NATIVE" ]; then python3 "$@";
        else docker run --rm -v "$PWD":/work -w /work pysecdec:1mass python3 "$@"; fi; }
IDS="B0_m1_0_m1 B0_m2_0_m2 B0_s_m1_m2 B0_Q_m1_m2 B0_t_0_m1 B0_u_0_m2 \
     C0_0_m1_t C0_0_m2_u C0_0_s_Q_a C0_0_s_Q_b C0_m1_m2_s C0_m1_u_Q C0_m2_t_Q D0_a D0_b"
for id in $IDS; do psd psd_fast_2mass.py build $id; done
echo "[build] all 15 master libraries built (kinematics = runtime real_parameters)"
