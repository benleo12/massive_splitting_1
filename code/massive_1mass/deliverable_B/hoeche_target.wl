(* Hoeche's MASSLESS spin remainder <P_{q->q}^{(1,p)}>, Eq.(2705) of 2505.10408,
   assembled from the explicit pieces (mu=1, s12=s, z the momentum fraction):
     c_Gamma = (4pi)^eps Gamma(1+eps)Gamma(1-eps)^2/Gamma(1-2eps)
     f2      = -c_Gamma/eps^2
     f1(z)   = -2 c_Gamma/eps^2 * 2F1(1,1;1-eps;1-z)          [Eq. line 2565]
     Pf      = C_F (1-eps)(1-z)        [fermionic tree, Eq.2606 with z2/z12=(1-z)]
     Pf1     = C_F (1-eps(1-z))        [new one-loop spin structure, line 2712]
   Remainder:
     R(z) = C_A (-1/s)^eps [ (1-z)f1(1-z) - (1/N^2)(z f1(z) - 2 f2) ] Pf
          - C_A (1+1/N^2) (-1/s)^eps  eps^2/(1-2eps) f2  Pf1 .
   This is the TARGET for the massive R^mag in the m->0 collinear limit. *)
NN=3; CF=(NN^2-1)/(2 NN); CA=NN;
cG = (4 Pi)^ep Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
f2 = -cG/ep^2;
f1[zz_] := -2 cG/ep^2 Hypergeometric2F1[1,1,1-ep,1-zz];
Pf[zz_]  := CF (1-ep)(1-zz);
Pf1[zz_] := CF (1-ep(1-zz));
brack[zz_] := (1-zz) f1[1-zz] - (1/NN^2)(zz f1[zz] - 2 f2);
Rem[zz_, sv_] := CA (-1/sv)^ep brack[zz] Pf[zz]
   - CA (1+1/NN^2)(-1/sv)^ep ep^2/(1-2ep) f2 Pf1[zz];
Print["Hoeche massless spin remainder <P_q^{(1,p)}>(z), as eps-series (s=1):"];
Do[Module[{r=Series[Rem[zv,1],{ep,0,2}]},
   Print["  z=",zv,":"];
   Do[Print["     eps^",k," : ", Chop[N[SeriesCoefficient[r,k],10],10^-9]],{k,-2,1}];
  ],{zv,{1/3, 1/2, 2/3}}];
