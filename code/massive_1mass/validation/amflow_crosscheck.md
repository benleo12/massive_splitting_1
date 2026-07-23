# AMFlow ↔ pySecDec cross-check (independent finite-mass oracle)

The 13 one-mass masters, computed two completely independent ways at a random
Euclidean point **s=−3, t=−5, u=−2, mq²=1** (Q²=−12):

- **pySecDec** — sector decomposition + quasi-Monte-Carlo (`run_masters.sh`).
- **AMFlow** — auxiliary-mass flow + differential equations + **Kira** IBP
  (`amflow_masters_1mass.wl`). Nothing shared with pySecDec.

Both in the same bare convention (ratio = 1, established on the bubble to 19 digits).
Reproduce: `wolframscript -file amflow_masters_1mass.wl` then
`wolframscript -file compare_amflow_psd.wl`.

## Result — max |AMFlow − pySecDec| per master (over ε⁻²…ε³)

| master | max\|Δ\| | master | max\|Δ\| |
|--------|--------|--------|--------|
| B0_s_mq_mq  | 9.3e−15 | C0_0_s_Q   | 3.7e−14 |
| B0_Q_mq_mq  | 9.8e−15 | C0_mq_u_Q  | 2.1e−11 |
| B0_u_0_mq   | 9.0e−13 | C0_mq_t_Q  | 8.3e−11 |
| B0_t_0_mq   | 1.5e−12 | C0_0_mq_t  | 2.0e−10 |
| B0_mq2_0_mq | 7.7e−12 | C0_0_mq_u  | 3.9e−10 |
|             |         | C0_mq_mq_s | 6.3e−10 |
| D0_a        | 1.3e−6  | D0_b       | 3.4e−7  |

**Worst over all 13 masters, all orders: 1.3×10⁻⁶** (a box at ε³).

The agreement is 10⁻¹⁰–10⁻¹⁶ through ε¹, degrading to ~10⁻⁶ only at ε³ for the
boxes. That residual is **pySecDec's QMC precision at ε³** (AMFlow is the 20-digit
side), i.e. it quantifies the pySecDec error rather than any disagreement — the
same "sector-decomposition precision degrades at higher ε" effect discussed for the
FONLL tables. Every master is confirmed by an independent method.

## AMFlow setup used (this machine)

AMFlow's shipped IBP backends need LiteRed, which is broken headless under MMA 14.2.
The working headless path: **AMFlow `Kira` reducer via a Docker bridge** —
`/Users/user/bin/kira-docker` translates AMFlow's kira-2.x CLI (`-p4`, `-s…`) to the
kira-3.1 image (`--parallel`, `--set_value`), with FireFly (no Fermat). Configured in
`AMFlow/ibp_interface/Kira/install.m` (`$KiraExecutable = /Users/user/bin/kira-docker`).
FiniteFlow was also built (`~/bassi/finiteflow`) and works, but its AMFlow path still
needs LiteRed/LiteIBP, so Kira is the route used here.
