#!/usr/bin/env bash
# Build the patched geval ("geval_big") used to evaluate the GPLs.
#
# WHY PATCHED: stock handyG ships geval.f90 with `character(len=20) line`, which
# silently truncates input lines longer than 20 chars.  The 1-mass box needs
# weight-5 GPLs in two letters, whose numeric input lines are far longer than
# that -> truncation -> wrong/zero results.  geval_big.f90 (in this folder) is
# identical except `character(len=500) line`.
#
# PREREQUISITE: a built handyG (https://gitlab.com/mule-tools/handyG).
#   git clone https://gitlab.com/mule-tools/handyG && cd handyG
#   ./configure && make            # produces libhandyg.a + the Fortran .mod files
#
# USAGE:
#   ./build_geval.sh /path/to/handyG        # the handyG source/build root
# then either put geval_big on your $PATH, or
#   export HANDYG_GEVAL=$(pwd)/geval_big
set -euo pipefail

HG="${1:?usage: build_geval.sh /path/to/handyG-build-root}"
HERE="$(cd "$(dirname "$0")" && pwd)"

# locate the handyG module include dir and static lib (layout varies by build)
INC="$(dirname "$(find "$HG" -name 'handyg.mod' -print -quit)")"
LIB="$(find "$HG" -name 'libhandyg.a' -print -quit)"
: "${INC:?could not find handyg.mod under $HG (build handyG first)}"
: "${LIB:?could not find libhandyg.a under $HG (build handyG first)}"

echo "include: $INC"
echo "lib:     $LIB"
gfortran -I"$INC" "$HERE/geval_big.f90" "$LIB" -o "$HERE/geval_big"
echo "built: $HERE/geval_big"
echo 'set it with:  export HANDYG_GEVAL='"$HERE"'/geval_big'
