# Deliverable B — massive magnetic (spin) remainder

The magnetic remainder completes the massive quark splitting function:

    R^mag  =  2Re<M0_q | M1_q> / |M0_q|^2  -  2Re<M0_sq | M1_sq> / |M0_sq|^2

i.e. the (quark − scalar) interference. It is the spin-dependent piece left over
once the semi-classical scalar dipole radiator (Deliverable A) is subtracted,
following the vertex decomposition of arXiv:2505.10408
(quark current = scalar + magnetic σ^{μν} + EOM).

## Status

All **ingredients are computed and validated**:
- the quark radiator (19 pure-QCD diagrams incl. n_f loops) reduced to 22 masters,
  all masters validated (AMFlow + analytic);
- the squark radiator (29 diagrams; the 6 closed-squark-loop diagrams shown to be
  separately gauge invariant and dropped to match the quark loop content);
- the reduced amplitudes are `../ampRaw_quark_1mass.txt`, `../ampRaw_squark_1mass.txt`.

What remains is the **final assembly + massless cross-check**, which is a bounded
computation, blocked only by local RAM (see below).

## The scheme and the physics

Following arXiv:2505.10408, the splitting functions are taken **unrenormalised**,
in a mixed gauge. The (quark − scalar) subtraction removes the semi-classical part;
the finite gauge dependence of the individual pieces is an external-leg
(A0(m^2), B0(m^2;0,m^2)) ambiguity that **vanishes as m→0**, so the massless limit
is gauge-clean.

Crucially, our amplitude regulates the q–g **collinear** region with the mass m_q
(this is the paper's central result: the mass screens the collinear singularity),
whereas the massless reference uses dimensional regularisation. Höche's massless
soft×collinear double pole therefore appears here as
`(1/eps_soft) × log(m_q^2)`, not as a pure `1/eps^2`.

## The target (locked)

Höche's explicit massless spin remainder ⟨P_{q→q}^{(1,p)}⟩ has

    ⟨P_{q→q}^{(1,p)}⟩ |_{eps^-2}  =  -2 C_A C_F (1-z)  =  -8(1-z) ,   z = s/(s+u) .

Verified: −4, −16/5, −8/3, −2 at z = 1/2, 3/5, 2/3, 3/4.
`hoeche_target.wl` / `target_table.wl` recompute this from Höche's f_1, f_2.

Collinear map: the q→qg splitting is the limit `t→0`, with `s_12 = t`, `z = s/(s+u)`.

## Running the cross-check

`rmag_massless_check.wl` assembles R^mag term-by-term (memory-frugal: each term is
spin-summed, then masters are substituted numerically and only the ε-coefficients
accumulate), at integer kinematics with the hierarchy `m_q^2 << |t| << |s|`
(`m_q^2=1` fixed, hard scales large) to expose the collinear `log(m_q^2/hard)`.
Two hard scales give the log slope, to be compared with `-8(1-z)`.

    wolframscript -file rmag_massless_check.wl

### Resource note

This step needs **~6 GB of free RAM**: the quark spin sum (`FermionSpinSum` over a
2.6M-leaf reduced amplitude) exceeds ~2 GB even term-by-term, because the Dirac
trace of a long γ-chain blows up intermediates. On a RAM-constrained machine it is
OOM-killed. Run it where memory is available (a workstation, or Perlmutter), or
supply the masters numerically via pySecDec at one collinear point and assemble
the interference once. The physics and the target above are settled; this is the
remaining mechanical evaluation.
