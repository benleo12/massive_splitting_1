# Capstone validation — massive reduction vs. validated massless amplitude

**Test:** `amplitudeMass1Ben(m²→0)` must reproduce `amplitudeMass0Stefan`
(the massless amplitude already validated digit-for-digit against Stefan's
notebook), order by order in ε. This validates the **freshly-generated FeynCalc
massive reduction** `ampRaw_1mass.txt` end-to-end (reduction + master-limit map +
normalization), independent of any per-master numeric oracle.

Reproduce: `wolframscript -file capstone_m0.wl`

## Result — ratio amp1(m→0)/amp0Stefan, per order ε⁻²…ε³

| (s,t,u)        | ε⁻² | ε⁻¹ | ε⁰  | ε¹  | ε²  | ε³  |
|----------------|-----|-----|-----|-----|-----|-----|
| (−1,−1,−1)     | 1   | 1   | 1   | 1.0000000 | 1.0000000 | 1   |
| (−2,−1,−1)     | 1   | 1   | 1   | 1   | 1   | 1   |
| (−1,−2,−½)     | 1   | 1   | 1   | 1   | 1   | 1   |
| (−3,−1,−½)     | 1   | 1   | 1   | 1   | 1   | 1   |
| (−1,−5,−2)     | 1   | 1   | 1   | 1   | 1   | 1   |

Three of the five points are **fully asymmetric** (s,t,u all distinct), where the
box masters contribute non-trivially (the box cancels only at s=u). Match to ε³ at
all of them ⇒ the reduction, the master→massless-limit map, **and the box arg
assignment** are all correct.

## Pinned facts

- **Box arg map** (the one piece that needed fixing):
  - `D0_a = D0(0,u,mq²,s; mq²,Q²; mq²,mq²,0,mq²)` → `boxF[u, s, Q²]`
  - `D0_b = D0(mq²,t,Q²,s; 0,mq²; mq²,0,mq²,mq²)` → `boxF[s, t, Q²]`
  with the **exact** `hypE` (₂F₁ expanded to ORD+3): the box is `1/ε²·[…]`, so the
  box's ε¹ needs the ε³ term of `hypE` — a truncated `1+ε²Li₂` fails from ε¹ at
  asymmetric points.
- Convention for the massless-limit forms: `cG / i / μ²` (identical to the package
  `boxF`/`StefanRules`), global `MSBarFac = (16π²/i)(e^γ/4π)^ε` and
  `cGRule` applied at assembly.
- Reduction coefficients are rational in (s,t,u,mq²) only — **no leftover D/ε** in
  the coefficients; all ε-dependence lives in the masters.
