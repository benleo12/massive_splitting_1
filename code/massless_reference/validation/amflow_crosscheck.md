# AMFlow ↔ pySecDec cross-check — two-mass masters (independent oracle)

The 15 two-mass masters, two independent ways at **s=−3, t=−5, u=−2, m1²=1, m2²=2**
(Q²=−13): **pySecDec** (sector decomposition + QMC) vs **AMFlow** (auxiliary-mass
flow + DEs + Kira IBP). Same bare convention. Reproduce:
`wolframscript -file amflow_masters_2mass.wl` then `... compare_amflow_psd_2mass.wl`.

## Result — max |AMFlow − pySecDec| per master (over ε⁻²…ε³)

| master | max\|Δ\| | master | max\|Δ\| |
|--------|--------|--------|--------|
| B0_s_m1_m2  | 8.9e−15 | C0_0_s_Q_b  | 3.0e−14 |
| B0_Q_m1_m2  | 9.3e−15 | C0_0_s_Q_a  | 3.5e−14 |
| B0_t_0_m1   | 1.5e−13 | C0_m2_t_Q   | 2.7e−11 |
| B0_u_0_m2   | 2.5e−13 | C0_m1_u_Q   | 3.1e−11 |
| B0_m2_0_m2  | 1.5e−11 | C0_m1_m2_s  | 1.9e−10 |
| B0_m1_0_m1  | 2.8e−11 | C0_0_m2_u   | 3.2e−10 |
|             |         | C0_0_m1_t   | 5.9e−10 |
| D0_a (2mass)| 6.1e−7  | D0_b (2mass)| 2.8e−7  |

**Worst over all 15 masters, all orders: 6.1×10⁻⁷** (a D0_2mass box at ε³).

The genuinely two-scale masters — the unequal-mass bubbles `B0(s;m1²,m2²)`,
`B0(Q²;m1²,m2²)`, the two orderings of `C0(0,s,Q²)`, the `C0(m1²,m2²,s)`, and the
**`D0_2mass`** boxes with two different internal masses — all agree with pySecDec to
10⁻¹⁰–10⁻¹⁵ through ε¹, limited only by pySecDec's QMC precision (~10⁻⁷) at ε³.
Independent confirmation of every two-mass master.

(AMFlow setup: same Docker-kira bridge as the 1-mass package —
see `../for_stefan_1mass/validation/amflow_crosscheck.md`.)
