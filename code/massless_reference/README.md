# 1-loop g\* → squark + gluon + antisquark, **two different masses** (numerical)

Third rung of the ladder, after `for_stefan/` (massless) and `for_stefan_1mass/`
(one mass). Here the two external squarks carry **different masses `m1 ≠ m2`**
(two heavy-quark flavours). Same `k1·ε(k2)` abelian projection, colour-stripped,
bare MS̄, delivered **numerically** (pySecDec), as an ε-Laurent series {ε⁻²…ε³}.

## The object

Two different external masses on a single squark line require a **flavour-off-diagonal
`S[1]` current** (`S̃₁–S̃₂*`): the model `feyncalc_reproduction/SQCDBGF_2mass`
has `MQU[1]=m1`, `MQU[2]=m2`, and the generation δ removed from the `S[1]` vertex.
The gluon vertices preserve flavour, so the squark line carries `m1` on the k1 side
of the current vertex and `m2` on the k3 side. Kinematic constraint:
**`s + t + u = Q² + m1² + m2²`**.

Same 6 amputated 1PI topologies as the 1-mass case; the fresh TID yields **15
masters** (vs 13 for one mass) — the extra ones are degeneracies split by the two
scales:

| kind | count | notable |
|------|-------|---------|
| bubbles | 6 | `B0(s; m1²,m2²)`, `B0(Q²; m1²,m2²)` — unequal-mass |
| triangles | 7 | `C0(0,s,Q²; m1²,m1²,m2²)` **and** `C0(0,s,Q²; m2²,m2²,m1²)` |
| boxes | 2 | `D0_2mass` — two different internal masses |

`m2 → m1` collapses the 15 masters **exactly** onto the 13 one-mass masters.

## What it computes

```
amplitudeMass2Ben[masterData, s, t, u, m1sq, m2sq]  ->  {ε⁻² .. ε³}
```
`masterData` = the 15 masters evaluated numerically at the point (pySecDec).
Reduction fixed in `ampRaw_2mass.txt`; assembler applies the loop measure `i` and
`MSBarFac` (identical convention to the 1-mass package).

### Quick start — FAST path (recommended)
Compile the masters **once** (kinematics are runtime parameters), then every point
is a ~2–3 s library call (no recompile):
```bash
cd validation
./build_masters_2mass.sh                                # ONE-TIME, ~2.5 min
./run_masters_fast_2mass.sh -2 -3 -1 1 2 masters.m      # ~2-3 s per point
cd ..
AMP_MASTERS=validation/masters.m wolframscript -file evaluate_amp_2mass.wl
```

### Quick start — simple path (recompiles each point, ~2.5 min/point)
```bash
cd validation && ./run_masters_2mass.sh -2 -3 -1 1 2 masters.m && cd ..
AMP_MASTERS=validation/masters.m wolframscript -file evaluate_amp_2mass.wl
```

No `AMP_MASTERS` → uses the bundled `validation/masters_demo_2mass.m` (m1²=1, m2²=2),
so the assembly runs without Docker.

**Timing (per point):** the cost is C++ codegen+compile, not integration (the QMC
itself is ~0.25 s/master). FAST path: ~2.5 min one-time build, then ~2–3 s/point;
in an in-process MC (no Docker restart) it's ~sub-second/point.

## Validation

1. **Capstone 1 — `m1,m2 → 0` reproduces the massless amplitude, to ε³.**
   `validation/capstone_m0_2mass.wl`: ratio = 1.000000 at every order ε⁻²…ε³ at
   5 points (3 fully asymmetric). Validates the 2-mass reduction end-to-end against
   the notebook-validated massless amplitude. (No pySecDec needed — closed massless
   master forms.)
2. **Structural nesting — `m2 → m1` → the 13 one-mass masters** (exact, from the
   generation).
3. **Capstone 2 — `m2 → m1` numerical:** `amplitudeMass2Ben(m1²=m2²=1)` ==
   `amplitudeMass1Ben(mq²=1)` at the same asymmetric point (box active).
   `validation/capstone_m1_2mass.wl`.
4. **Pole structure:** finite with `ε⁻² = 0` (massive legs screen the collinear
   singularity; only the soft-gluon ε⁻¹ survives) — at `m1≠m2` too.
5. **Independent AMFlow cross-check:** all 15 masters (incl. the `D0_2mass` boxes)
   confirmed by AMFlow (aux-mass flow + Kira IBP) vs pySecDec to 10⁻¹⁰–10⁻¹⁵ through
   ε¹, worst 6.1×10⁻⁷ at ε³ (pySecDec QMC). See `validation/amflow_crosscheck.md`.

### Reference numbers
```
amplitudeMass2Ben[-2,-3,-1, m1²=1, m2²=2]  (Q²=-9):
  ε⁻²: 0   ε⁻¹: 0.6401   ε⁰: -1.0291   ε¹: 2.5835   ε²: -1.7455   ε³: 5.2559
```

## Files
```
ampRaw_2mass.txt               frozen FeynCalc reduction (rational in s,t,u,m1²,m2²)
master_list_2mass.txt          the 15 masters
amplitude_functions_2mass.wl   assembler -> amplitudeMass2Ben[masterData, s,t,u,m1sq,m2sq]
evaluate_amp_2mass.wl          CLI driver
validation/
  psd_master_2mass.py          one master via pySecDec (two scales)
  run_masters_2mass.sh         all 15 masters at a point -> masters.m
  parse_psd.py                 pySecDec output -> Mathematica master data
  capstone_m0_2mass.wl         m1,m2->0 capstone (run it)
  capstone_m1_2mass.wl         m2->m1 numerical capstone
  capstone_results.md          capstone tables
  masters_demo_2mass.m         bundled masters for m1²=1, m2²=2
feyncalc_reproduction/         SQCDBGF_2mass model + 04mm generation script
                               (under for_stefan_1mass/feyncalc_reproduction/)
```

## Scope
Euclidean points (`s,t,u<0`, `m1²,m2²>0`) keep masters real. Out of scope, as
before: massive-quark **loop** insertions, real external quarks, non-abelian colour,
other Lorentz projections.
