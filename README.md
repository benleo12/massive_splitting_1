# Massive one-loop scalar dipole radiator

One-loop QCD splitting amplitudes for gluon emission off a pair of **massive**
colour charges, formulated as scalar dipole radiators — the massive counterpart
of the box-type contributions (Fig. 14) of

> J. M. Campbell, S. Höche, M. Knobbe, C. T. Preusz, D. Reichelt,
> *QCD splitting functions beyond kinematical limits*, arXiv:2505.10408.

Everything here is **fully numerical**: each master integral is generated once by
sector decomposition with the invariants and masses kept as runtime parameters,
so an arbitrary phase-space point costs a fraction of a second. The two heavy-quark
masses are kept independent (equal-mass and unequal-mass dipoles both covered).

---

## Repository layout

| directory | contents |
|---|---|
| `paper/` | The write-up (JHEP style). `main.tex` + appendix coefficient files `app_*.tex`. |
| `amplitudes/` | The deliverable **reduced amplitudes** as plain text, in terms of Passarino–Veltman masters (colour/coupling stripped). Abelian (`ampRaw_1mass`), non-abelian (`ampRaw_nab_1mass`), squark (`ampRaw_squark_1mass`); `master_list_1mass.txt` names the basis. |
| `code/massive_1mass/` | The massive one-mass calculation: FeynArts/FeynCalc generation + reduction (`feyncalc_reproduction/`), the amplitude assemblers (`amplitude_functions_*.wl`), the master closed-forms and pySecDec builders, and all validation scripts (`validation/`). |
| `code/massless_reference/` | The massless reference calculation (Höche's massless result reproduced), used to validate the `m→0` limit. |
| `code/massive_2mass/` | The unequal-mass ($m_i\neq m_k$) calculation: model with independent masses per generation, 15-master basis incl. the genuine two-mass box, validated by $m_k\to m_i$ collapse onto the one-mass amplitude. |
| `collaborator_notebook/` | S. Höche's original massless one-loop notebook `one-loop_35.nb` and the `SQCDBGF` FeynArts model it uses (standalone). |

pySecDec build directories (`m_*/`, the compiled `.so` libraries) are **not**
committed — they are platform-specific and regenerable; see *Reproducing* below.

---

## What is verified

The scalar dipole radiator (Deliverable A) is validated three independent ways.

1. **Every master integral vs auxiliary-mass flow (AMFlow).**
   All 22 masters of the equal-mass radiator agree with an independent AMFlow
   evaluation at a random Euclidean point to the quasi-Monte-Carlo precision
   (`code/massive_1mass/validation/amflow_*`, `capstone_results.md`).
   Bubbles/triangles agree to 10⁻⁹–10⁻¹⁶; the boxes to ~10⁻⁷ at ε³.

2. **Massless limit vs the closed-form massless expressions.**
   - *Abelian:* the assembled amplitude sent to `m→0` reproduces the massless
     result order by order through ε³ (ratio = 1 to 10 digits, several points
     incl. fully asymmetric — `validation/capstone_m0.wl`).
   - *Non-abelian:* dividing the `m→0` limit by the massless non-abelian
     combination gives a **kinematics-independent** ratio `−16π²i (1 + 0.5675 ε + …)`,
     identical at four points to ~18 digits (`feyncalc_reproduction/07nab_verify_m0.wl`).

3. **Soft limit vs the one-loop massive soft current** of Bierenbaum–Czakon–Mitov
   (arXiv:1107.4384, Case 3). The double pole is reproduced with a single
   kinematics-independent normalisation, `−64π²i`, at eight physical phase-space
   configurations (two dipole invariants × a decade in softness), 11–13 digits
   (`validation/bcm_soft_compare.wl`).

The radiator additionally satisfies the QCD Ward identity **symbolically**
(the reduced amplitude carries the exact dipole structure
`p_i/(t−m²) − p_k/(u−m²)`; verified as an algebraic identity).

---

## Status of the two deliverables

- **A — massive scalar dipole radiator + the no-sectorization argument.**
  Complete and validated (above). Written up in `paper/`. The physics point:
  the quark mass screens the collinear singularity, so the massive radiator has
  no abelian `1/ε²` double pole and its entire singular structure is soft;
  consequently the massive NNLO subtraction needs no phase-space sectorization.

- **B — massive magnetic (spin) remainder** = (quark − scalar) interference.
  *Not needed for the heavy-quark subtraction* (Czakon): the remainder is singular
  only in the collinear limit, which the mass regulates, so it is finite and
  subtraction-irrelevant; the scalar radiator is the complete ingredient. The
  ingredients (quark radiator, 19 diagrams/22 masters; squark radiator, 29
  diagrams — all masters validated) and the massless cross-check target
  `⟨P_{q→q}^{(1,p)}⟩|_{ε⁻²} = −2 C_A C_F (1−z)`, `z=s/(s+u)`, are kept in
  `code/massive_1mass/deliverable_B/` for reference.

---

## Reproducing

Requirements: Mathematica with **FeynArts + FeynCalc 9.3.1** (10.1.0 aborts on
this model), **pySecDec** (native `pip install pySecDec`, or the Docker image),
and optionally **AMFlow** for the independent master cross-check.

```bash
cd code/massive_1mass

# 1. build the master-integral libraries once (kinematics as runtime parameters)
bash validation/build_masters.sh          # ~minutes; creates validation/m_*/

# 2. evaluate any phase-space point (sub-second per point thereafter)
python3 validation/psd_fast.py evalall  <s> <t> <u> <mq2>

# 3. assemble the amplitude (masters -> Laurent series in eps)
wolframscript -file amplitude_functions_nab.wl     # non-abelian radiator
wolframscript -file amplitude_functions_1mass.wl   # abelian radiator

# 4. validation
wolframscript -file validation/capstone_m0.wl      # m->0 vs massless (ratio = 1)
wolframscript -file validation/bcm_soft_compare.wl # soft limit vs BCM (-64 pi^2 i)
```

The reduced amplitudes in `amplitudes/` are the direct output of the FeynCalc
reduction (`feyncalc_reproduction/`), colour- and coupling-stripped; the explicit
reduction coefficients are also tabulated in the paper appendix (`paper/app_*.tex`).

## Conventions

- Kinematics: `s=(p_i+p_k)²`, `t=(p_i+q)²`, `u=(p_k+q)²`, `Q²=(p_i+p_k+q)²`,
  with `p_i²=p_k²=m²`, `q²=0`, so `s+t+u = Q²+2m²`.
- Master normalisation and the `c_Γ`/loop-measure factors follow Ref. 2505.10408;
  see the header of `amplitude_functions_1mass.wl` and the paper.
