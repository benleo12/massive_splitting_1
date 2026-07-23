# pySecDec reference values for the 13 massive masters

Point: s=t=u=−1, mq²=1, Q²=s+t+u−2mq²=−5. pySecDec (Docker, `pysecdec:1mass`),
ε⁻¹…ε³, default Γ(ε) loop prefactor (same convention as the 2F1 bubble forms).
This both (a) validates the Gate-1 master basis and (b) is the target for the
symbolic DE work. Symmetry degeneracies at t=u confirm the reduction.

## Bubbles (DONE — symbolic 2F1, matched to 16 digits; see masters_bubbles.wl)
- B0(t,0,mq²)=B0(u,0,mq²):  ε⁻¹ 1, ε⁰ 0.03649, ε¹ 1.17820, ε² −0.09516, ε³ 1.17148  (at p²=−1)
- B0(mq²,0,mq²):            (p²=mq²=1 case)
- B0(s,mq²,mq²)=B0(Q²,mq²,mq²): equal-mass 2F1

## Triangles (pySecDec refs; symbolic TODO)
- C0_0_mq_t = C0_0_mq_u (finite): ε⁰ −1.2337006, ε¹ −1.3914883, ε² −4.0646804, ε³ −6.6544582
- C0_0_s_Q  (finite, all-massive-internal): ε⁰ −0.3473472, ε¹ 0.3200055, ε² −0.4400026, ε³ 0.4541364
- C0_mq_mq_s (1/ε pole): ε⁻¹ 0.4304089, ε⁰ −0.3119964, ε¹ 0.4680400, ε² −0.4570732, ε³ 0.4860404
- C0_mq_t_Q = C0_mq_u_Q (finite): ε⁰ −0.9293702, ε¹ −1.0502411, ε² −3.6010698, ε³ −6.1765686

## Boxes (pySecDec refs; symbolic TODO — the crux)
- D0_a = D0_b (1/ε pole, equal at t=u): ε⁻¹ −0.2152045, ε⁰ −0.2312493, ε¹ −0.6617568, ε² −1.4394271, ε³ −3.1886355

Reproduce any: `docker run --rm -v "$PWD":/work pysecdec:1mass python3 psd_master.py <id> <s> <t> <u> <mq2>`
ids: B0_t_0_mq B0_u_0_mq B0_mq2_0_mq B0_s_mq_mq B0_Q_mq_mq
     C0_0_mq_t C0_0_mq_u C0_0_s_Q C0_mq_mq_s C0_mq_t_Q C0_mq_u_Q D0_a D0_b
