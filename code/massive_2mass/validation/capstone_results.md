# 2-mass capstone validations

The two external squarks carry different masses `m1 ≠ m2` (flavor-off-diagonal
`S[1]` current; gluon vertices preserve flavor, so the squark line carries `m1` on
the k1 side of the current vertex and `m2` on the k3 side). Constraint
`s+t+u = Q² + m1² + m2²`. 15 masters (6 B0, 7 C0, 2 D0), two-scale, including the
genuine `D0_2mass` boxes (two different internal masses).

## Capstone 1 — `m1,m2 → 0` reproduces the massless amplitude (to ε³)

`capstone_m0_2mass.wl`: send both masses to zero, every 2-mass PaVe → its massless
limit (`Bf/Cf1/Cf2/boxF`), compare to `amplitudeMass0Stefan` (notebook-validated).

| (s,t,u)     | ε⁻² | ε⁻¹ | ε⁰ | ε¹ | ε² | ε³ |
|-------------|-----|-----|----|----|----|----|
| (−1,−1,−1)  | 1   | 1   | 1  | 1.0000000 | 1.0000000 | 1 |
| (−2,−1,−1)  | 1   | 1   | 1  | 1  | 1  | 1 |
| (−1,−2,−½)  | 1   | 1   | 1  | 1  | 1  | 1 |
| (−3,−1,−½)  | 1   | 1   | 1  | 1  | 1  | 1 |
| (−1,−5,−2)  | 1   | 1   | 1  | 1  | 1  | 1 |

Ratio = 1.000000 at every order, 5 points (3 fully asymmetric). Validates the
2-mass reduction `ampRaw_2mass.txt` end-to-end (reduction + master-limit map +
colour strip + box assignment).

## Structural nesting — `m2 → m1` collapses the master set to the 13 one-mass masters

Set `m2=m1` in the master argument list: the 15 two-mass masters collapse exactly
to the 13 one-mass masters (the two `C0(0,s,Q²)` orderings merge, the two boxes →
the one-mass boxes, etc.). Confirms the 2-mass amplitude nests onto the validated
1-mass one. (See 04mm generation log.)

## Capstone 2 — `m2 → m1` numerical (amplitudeMass2Ben == amplitudeMass1Ben)

`capstone_m1_2mass.wl`: pySecDec the 15 masters at `m1²=m2²=1`, assemble
`amplitudeMass2Ben`, compare to `amplitudeMass1Ben(mq²=1)` at the SAME asymmetric
point `s=-2,t=-3,u=-1` (Q²=-8, box active).

```
amplitudeMass2Ben(m1=m2=1) : {0, 0.673564, -0.743165, 2.138792, -0.441816, 4.533664}
amplitudeMass1Ben(mq2=1)   : {0, 0.673564, -0.743165, 2.138792, -0.441816, 4.533669}
|diff| eps^-2..3           : {0, 2e-14, 1e-9, 6e-8, 4e-8, 5e-6}
```
Agreement to pySecDec QMC precision (the two sides are independent pySecDec runs:
15 two-mass masters vs 13 one-mass masters), degrading gently from ε⁻¹ (2e-14) to
ε³ (5e-6). Confirms the 2-mass reduction + 15-master numerics + assembly nest onto
the validated 1-mass pipeline.
