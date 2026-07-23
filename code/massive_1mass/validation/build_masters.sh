#!/bin/bash
# ONE-TIME build of all 13 master libraries with kinematics as runtime parameters.
# Uses native pySecDec if importable (pip install pySecDec), else the Docker image.
set -e
cd "$(dirname "$0")"
if python3 -c "import pySecDec" 2>/dev/null; then NATIVE=1; fi
psd() { if [ -n "$NATIVE" ]; then python3 "$@";
        else docker run --rm -v "$PWD":/work -w /work pysecdec:1mass python3 "$@"; fi; }
IDS="B0_t_0_mq B0_u_0_mq B0_mq2_0_mq B0_s_mq_mq B0_Q_mq_mq \
     C0_0_mq_t C0_0_mq_u C0_0_s_Q C0_mq_mq_s C0_mq_t_Q C0_mq_u_Q D0_a D0_b"
for id in $IDS; do psd psd_fast.py build $id; done
echo "[build] all 13 master libraries built (kinematics = runtime real_parameters)"
