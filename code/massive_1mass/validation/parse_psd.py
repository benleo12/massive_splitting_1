#!/usr/bin/env python3
# Parse PSDRESULT lines (pySecDec with-prefactor output) from one or more files
# into a Mathematica association  <| "id" -> {kmin, {re_kmin, ..., re_3}} |>
# central values only (imag parts are ~1e-14 numerical noise at Euclidean points).
# Usage:  python3 parse_psd.py out.m file1.txt [file2.txt ...]
import sys, re

out_m = sys.argv[1]
files = sys.argv[2:]
term = re.compile(r'\(\(([-+0-9.eE]+),[^)]*\)\s*\+/-\s*\([^)]*\)\)(\*eps(?:\^(-?\d+))?)?')

pt = re.compile(r's=([-0-9./]+)\s+t=([-0-9./]+)\s+u=([-0-9./]+)\s+mq2=([-0-9./]+)')
masters = {}
point = None
for fn in files:
    for line in open(fn):
        if not line.startswith('PSDRESULT'):
            continue
        if point is None:
            pm = pt.search(line)
            if pm:
                point = pm.groups()
        mid = line.split()[1]
        body = line.split('::', 1)[1]
        coeffs = {}
        for m in term.finditer(body):
            re_val = float(m.group(1))
            if m.group(2) is None:            # no *eps  -> eps^0
                p = 0
            elif m.group(3) is None:          # *eps     -> eps^1
                p = 1
            else:
                p = int(m.group(3))           # *eps^N
            coeffs[p] = re_val
        if coeffs:
            masters[mid] = coeffs

def mma_num(x):
    return repr(x).replace('e', '*^')

with open(out_m, 'w') as f:
    f.write('<|\n')
    if point:
        f.write('  "POINT" -> {%s},\n' % ', '.join(point))   # {s, t, u, mq2}
    for i, (mid, c) in enumerate(sorted(masters.items())):
        kmin = min(c)
        lst = [c.get(k, 0.0) for k in range(kmin, 4)]
        comma = ',' if i < len(masters) - 1 else ''
        f.write('  "%s" -> {%d, {%s}}%s\n' % (mid, kmin, ', '.join(mma_num(v) for v in lst), comma))
    f.write('|>\n')
print("wrote %s : %d masters %s" % (out_m, len(masters), sorted(masters)))
