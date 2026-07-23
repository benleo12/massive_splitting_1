# 1-loop g\* → squark + gluon + antisquark, **massive squark line** (numerical)

Companion to `for_stefan/` (the massless cross-check). Same object — the off-shell
scalar → squark + background-gluon + antisquark amplitude `S[1] → S[13] V[50] -S[13]`,
abelian (1/N) colour, `k1·ε(k2)` projection, colour/structure-stripped — but with the
**squark line made massive** (`MQU = m`, external squark shells `k1²=k3²=m²`, gluon
`k2²=0`, so `s+t+u = Q² + 2m²`). Same 6 amputated 1PI diagrams, bare MS̄, no
wavefunction/mass counterterms.

Stefan never computed this analytically, so it is delivered **numerically** — in line
with the all-numerical (pySecDec / Czakon-pragmatic / STRIPPER-style) strategy: the
masses are just numbers in a sector-decomposition integral, and the same evaluator
plugs into a fully-differential MC as the virtual building block. Higher ε orders are
numerical and the precision is set per point by the QMC.

## What it computes

```
amplitudeMass1Ben[masterData, s, t, u, mq2]  ->  {ε⁻² .. ε³}
```

`masterData` is the set of 13 master integrals evaluated **numerically at this point**
(pySecDec). The reduction (rational coefficients in `s,t,u,mq²`) is fixed and lives in
`ampRaw_1mass.txt`; the assembler multiplies in the masters, the loop measure factor
`i`, and the global `MSBarFac = (16π²/i)(e^γ/4π)^ε`, and returns the ε-Laurent series.

### Quick start — FAST path (recommended)
Compile the 13 masters **once** (kinematics are runtime parameters), then every
point is a ~2–3 s library call:
```bash
cd validation
./build_masters.sh                              # ONE-TIME, ~2 min
./run_masters_fast.sh -2 -3 -1 1 masters.m      # ~2-3 s per point
cd ..
AMP_MASTERS=validation/masters.m wolframscript -file evaluate_amp_1mass.wl
```

### Quick start — simple path (recompiles each point, ~2 min/point)
```bash
cd validation && ./run_masters.sh -2 -3 -1 1 masters.m && cd ..
AMP_MASTERS=validation/masters.m wolframscript -file evaluate_amp_1mass.wl
```

With no `AMP_MASTERS` set, the driver uses the bundled `validation/masters_demo_asym.m`,
so the assembly runs without Docker. Timing is dominated by C++ compile, not the QMC
integration (~0.25 s/master) — hence the one-time build + fast per-point path.

## Master integrals (13)

| kind | ids | how |
|------|-----|-----|
| bubbles (5) | `B0_{t,u}_0_mq`, `B0_mq2_0_mq`, `B0_{s,Q}_mq_mq` | exact 2F1 / closed form (`masters_bubbles.wl`) **and** pySecDec |
| triangles (6) | `C0_0_mq_{t,u}`, `C0_0_s_Q`, `C0_mq_mq_s`, `C0_mq_{t,u}_Q` | pySecDec |
| boxes (2) | `D0_a`, `D0_b` | pySecDec |

Convention: each master is the **bare** integral `Γ(n−D/2) ∫ (Symanzik)` (identical to
pySecDec); the amplitude convention is reached by the universal loop-measure factor `i`
(proved analytically for the bubble, `Bf = i·B0zm`; confirmed by the pole structure for
the triangles/boxes — see below). `B0_mq2_0_mq` is the squark on its own shell `p²=m²`,
where the 2F1 form has a threshold pole; it has the clean closed form
`Γ(ε)(m²)^{−ε}/(1−2ε)` (`ε⁰ = 2−γ`, matches pySecDec).

## Validation (this is the evidence that replaces the absent analytic reference)

1. **Capstone — `m²→0` reproduces the validated massless amplitude, to ε³.**
   `validation/capstone_m0.wl` takes the freshly-generated massive reduction, sends
   `m²→0` (every massive PaVe → its massless limit), and compares to
   `amplitudeMass0Stefan` (already validated digit-for-digit vs Stefan's notebook).
   **Ratio = 1.000000 at every order ε⁻²…ε³ at 5 points (2 symmetric + 3 fully
   asymmetric).** This validates the reduction, the master→limit map, and the box arg
   assignment end-to-end. See `validation/capstone_results.md`.

2. **Masters two ways.** The 5 bubbles agree between the independent 2F1 closed forms
   and pySecDec to **~10⁻¹³** (pySecDec's QMC precision), including `B0_mq2_0_mq`.

3. **Pole structure / box convention.** The assembled amplitude is finite with
   **`ε⁻² = 0`** (the massive squark legs screen the collinear singularity; only the
   soft-gluon `ε⁻¹` pole survives) at both a symmetric and a fully asymmetric point.
   The asymmetric point (box contributing non-trivially) still cancels `ε⁻²` exactly —
   a wrong box normalization would spoil the cancellation, so this confirms the box
   convention.

4. **Independent AMFlow cross-check (finite mass).** All 13 masters, computed a second,
   fully independent way — **AMFlow** (auxiliary-mass flow + DEs + Kira IBP) vs pySecDec
   — agree to 10⁻¹⁰–10⁻¹⁶ through ε¹, limited only by pySecDec's QMC at ε³ (worst
   1.3×10⁻⁶). See `validation/amflow_crosscheck.md`. This is the independent finite-mass
   oracle (nothing shared with pySecDec's sector decomposition).

### Reference numbers

```
amplitudeMass1Ben[-1,-1,-1, mq2=1]  (Q2=-5, box cancels at s=u):
  ε⁻²: 0   ε⁻¹: 1.1942   ε⁰: -1.7655   ε¹: 4.2992   ε²: -2.5652   ε³: 8.5400

amplitudeMass1Ben[-2,-3,-1, mq2=1]  (Q2=-8, fully asymmetric, box active):
  ε⁻²: 0   ε⁻¹: 0.6736   ε⁰: -0.7432   ε¹: 2.1388   ε²: -0.4418   ε³: 4.5337
```

## Files

```
ampRaw_1mass.txt              frozen FeynCalc reduction (rational in s,t,u,mq²; PaVe B0/C0/D0)
master_list_1mass.txt         the 13 masters
amplitude_functions_1mass.wl  the assembler -> amplitudeMass1Ben[masterData, s,t,u,mq2]
masters_bubbles.wl            exact 2F1 bubble masters (validation cross-check)
evaluate_amp_1mass.wl         CLI driver
validation/
  psd_master.py               one master via pySecDec
  run_masters.sh              all 13 masters at a point -> masters.m
  parse_psd.py                pySecDec output -> Mathematica master data
  capstone_m0.wl              the m²→0 capstone (run it)
  capstone_results.md         capstone table
  master_refs.md              pySecDec reference values (symmetric point)
  masters_demo_asym.m         bundled masters for s=-2,t=-3,u=-1,mq²=1
  Dockerfile                  builds the pysecdec:1mass image
feyncalc_reproduction/        how ampRaw_1mass.txt was generated (massive amps + TID)
```

## Notes / scope
- Euclidean points (`s,t,u<0, mq²>0`) keep all masters real and the QMC well-behaved.
- For a phenomenology MC the masters are evaluated at the points the MC visits; no
  re-derivation per mass/kinematics — the reduction and the assembler are fixed.
- Out of scope (as in the massless package): massive-quark **loop** contribution, real
  external quarks, non-abelian colour, other Lorentz projections.
