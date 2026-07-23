Independent FeynCalc reproduction of Stefan's amplitude reduction
=================================================================

**Top-to-bottom** verification: we regenerate Stefan's 1-loop
g* → q q̄ g amplitude reduction ourselves in FeynCalc — diagram
generation, Dirac + colour algebra, tensor reduction (`TID`) and
Passarino–Veltman reduction (`PaVeReduce`) — and confirm our reduced
amplitude equals his, *exactly*, for **both** colour projections.

This matters because for the massive case there is no reference; the
whole chain (FeynCalc → reduction → masters → evaluation) has to be
trusted, validated here on the massless case where Stefan's result exists.

The two colour structures
-------------------------

The amplitude splits into a **non-abelian** part (colour factor
`f^{abc}`, process `S[1]→S[13] V[50] -S[13]`, diagrams {14},{13,15},{5,12})
and an **abelian** part (the `1/N` structures, process
`S[1]→S[13] V[50] S[13]` = `diagsl2`, diagrams {1,2},{9,11},{5,8}).
Both are reproduced and verified.

| step | script | what it does |
|---|---|---|
| 1 | `01_nonabelian_generate_TID.wl`   | generate diagrams, colour projections, momenta routing, `TID` → `/tmp/nab_amp1.mx` |
| 2 | `02_nonabelian_reduce_assemble.wl`| `PaVeReduce` → scalar masters, divide by `nabprefac`, extract k1·ε coeff → `/tmp/nab_coeffK1.mx` |
| 3 | `03_nonabelian_verify.wl`         | compare to Stefan's hand-derived `MyBoxNAbA` (cell 40) via Package-X |
| 4 | `04_abelian_generate_reduce.wl`   | `diagsl2` + abelian colour projections (`TaTAbaRule`, 1/N rules), `TID`, `PaVeReduce`, assemble → `/tmp/ab_sum.mx` |
| 5 | `05_abelian_extract_k1.wl`        | extract the `SPD[k1,ε(k2)]` coefficient → `/tmp/ab_coeffK1.mx` |
| 6 | `06_abelian_verify.wl`            | compare to Stefan's hand-derived `MyBoxAmpAbLS1` (cell 90) via Package-X |

Result
------

For **both** colour structures, our FeynCalc-reduced amplitude equals
Stefan's hand-derived master combination to all orders in ε, at every
kinematic point:

    our_reduced_amplitude  =  N(ε) · Stefan_reference

with the **same** kinematic-independent, colour-structure-independent
normalisation (the Package-X ↔ LoopFuncRules convention factor):

    N(ε) = −16π² i + 89.62 i ε − 25.43 i ε² + …      (leading −16π² i = MSBarFac)

The residual `our − N(ε)·Stefan_reference` is **identically zero**
(exact symbolic — Package-X gives exact results) at ε⁻², ε⁻¹, ε⁰ for
(s,t,u) = (−1,−1,−1), (−2,−1,−1), (−1,−2,−½), (−3,−2,−2).

That N(ε) is identical for the abelian and non-abelian pieces is itself
a consistency check — it is the pure loop normalisation, the same for
every colour structure.

So our independent reduction reproduces Stefan's amplitude exactly.  The
*master integrals* in that reduction are then evaluated with our
HyperForm computation (`../evaluate_amp.wl`), which reproduces the
literature master values to 12+ digits.  Both halves of the pipeline —
the FeynCalc front-end and the HyperForm master evaluation — are
independently validated, on both colour projections.

For the massive case the same front-end is rerun with massive
propagators, and the HyperForm side with massive masters; every link in
the chain has been checked here where Stefan's answer exists.

Notes
-----

* FeynCalc 9.3.1 (the SQCDBGF model does not initialise under FeynArts
  shipped with FeynCalc 10.x).  Step 1/4 symlink the `FeynCalc backup`
  directory as `/tmp/fc93/FeynCalc`.
* FeynHelpers + Package-X needed for `PaXEvaluate` (steps 3, 6).
* Run order per colour structure: 1→2→3 (non-abelian), 4→5→6 (abelian);
  each caches to `/tmp`.
