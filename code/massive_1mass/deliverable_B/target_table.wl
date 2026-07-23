(* The massless cross-check TARGET, tabulated at z = s/(s+u) for kinematic points
   we can drive collinear (t->0).  This is what our massive R^mag must reproduce as
   m->0, t->0, with z=s/(s+u), s12=t. *)
NN=3; CF=(NN^2-1)/(2 NN); CA=NN;
cG=(4 Pi)^ep Gamma[1+ep]Gamma[1-ep]^2/Gamma[1-2 ep];
f2=-cG/ep^2; f1[zz_]:=-2 cG/ep^2 Hypergeometric2F1[1,1,1-ep,1-zz];
Pf[zz_]:=CF(1-ep)(1-zz); Pf1[zz_]:=CF(1-ep(1-zz));
brack[zz_]:=(1-zz)f1[1-zz]-(1/NN^2)(zz f1[zz]-2 f2);
Rem[zz_,s12_]:=CA(-1/s12)^ep brack[zz]Pf[zz]-CA(1+1/NN^2)(-1/s12)^ep ep^2/(1-2ep)f2 Pf1[zz];
(* pull out the universal (-1/s12)^eps and give the z-dependent coefficient series *)
Print["Hoeche remainder / (-1/s12)^eps , as eps-series (double pole is s12-clean):"];
Do[Module[{r=Series[Rem[zv,1]/(-1)^ep,{ep,0,1}]},   (* s12=1, strip the log later *)
  Print["  z=",N[zv,4]," (=s/(s+u)):  ",
    "eps^-2=",Chop[N[SeriesCoefficient[r,-2],8]],
    "  eps^-1=",Chop[N[SeriesCoefficient[r,-1],8]],
    "  eps^0=",Chop[N[SeriesCoefficient[r,0],8]]]],
 {zv,{1/2, 3/5, 2/3, 3/4}}];
Print[""];
Print["  Structure: eps^-2 coefficient vs z (the cleanest test target):"];
Do[Print["    z=",N[zv,4]," -> ", Simplify[SeriesCoefficient[Series[Rem[zv,1],{ep,0,-2}],-2]]],
 {zv,{1/2,3/5,2/3,3/4}}];
