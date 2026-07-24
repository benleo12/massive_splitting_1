# Massive factorizable set: Ward check + hard-vertex function

`ward_factorizable_massive.wl` computes, with the squark mass kept:
- the squark self-energy  Sigma(p^2) = (CF g^2/16pi^2)(A0(m^2) - 2 pbar^2 B0(p^2;0,m^2)),
  pbar^2 the inverse-propagator virtuality (reduces to the massless -2p^2 I2(p^2) as m->0);
- the squark-squark-gluon vertex correction with both scalar legs off shell,
  contracted with the gluon momentum.

Result: k.Gamma is proportional to Sigma(q^2)-Sigma(p^2) with a kinematics-independent
constant, verified numerically at three off-shell points. The massive factorizable
(Fig. 2 type) set therefore satisfies the naive background-field Ward identities,
so its abelian part cancels the self-energies exactly as in the massless case.

`ihf_massive.wl` computes the massive factorizable hard-vertex function (topology (e)
of the box figure), from the production vertex with the emitting leg off shell:

  I_hf^m(Q^2,t) = g^2/(16pi^2) Ti.Tk [ (3m^2-2Q^2+t) C0(m^2,Q^2,t;0,m^2,m^2)
                   - B0(Q^2;m^2,m^2) + B0(t;0,m^2) + B0(m^2;0,m^2) ]

whose m->0 limit reproduces the massless I_{h,f} of arXiv:2505.10408 exactly.
All four masters are in the validated basis (compiled libraries evaluate them).
