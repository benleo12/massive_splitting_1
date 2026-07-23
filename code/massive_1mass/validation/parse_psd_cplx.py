#!/usr/bin/env python3
# Like parse_psd.py, but KEEPS THE IMAGINARY PARTS.
# Needed for PHYSICAL (Minkowski) kinematics, where the masters genuinely have
# absorptive parts -- parse_psd.py discards them (fine only at Euclidean points).
# Usage:  python3 parse_psd_cplx.py out.m file1.txt [file2.txt ...]
import sys, re

out_m = sys.argv[1]
files = sys.argv[2:]
term = re.compile(r'\(\(([-+0-9.eE]+),([-+0-9.eE]+)\)\s*\+/-\s*'
                  r'\(([-+0-9.eE]+),([-+0-9.eE]+)\)\)(\*eps(?:\^(-?\d+))?)?')
pt = re.compile(r's=([-0-9.eE+/]+)\s+t=([-0-9.eE+/]+)\s+u=([-0-9.eE+/]+)\s+mq2=([-0-9.eE+/]+)')

masters, errs, point = {}, {}, None
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
        coeffs, ec = {}, {}
        for m in term.finditer(body):
            re_v, im_v, re_e, im_e = (float(m.group(i)) for i in (1, 2, 3, 4))
            if m.group(5) is None:      p = 0     # no *eps  -> eps^0
            elif m.group(6) is None:    p = 1     # *eps     -> eps^1
            else:                       p = int(m.group(6))
            coeffs[p] = (re_v, im_v); ec[p] = max(abs(re_e), abs(im_e))
        if coeffs:
            masters[mid] = coeffs; errs[mid] = ec

def mma(x):
    return repr(x).replace('e', '*^')

with open(out_m, 'w') as f:
    f.write('<|\n')
    if point:
        f.write('  "POINT" -> {%s},\n' % ', '.join(point))       # {s, t, u, mq2}
    items = sorted(masters.items())
    for i, (mid, c) in enumerate(items):
        kmin = min(c)
        lst = ['(%s + I*(%s))' % (mma(c.get(k, (0.0, 0.0))[0]), mma(c.get(k, (0.0, 0.0))[1]))
               for k in range(kmin, 4)]
        comma = ',' if i < len(items) - 1 else ''
        f.write('  "%s" -> {%d, {%s}}%s\n' % (mid, kmin, ', '.join(lst), comma))
    f.write('|>\n')

worst = max((max(e.values()), m) for m, e in errs.items()) if errs else (0, '-')
print("wrote %s : %d masters, worst abs MC error %.2e (%s)" % (out_m, len(masters), worst[0], worst[1]))
