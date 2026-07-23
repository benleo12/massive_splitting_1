(* ============================================================================
   CAPSTONE 1 for the 2-mass amplitude:  amplitudeMass2Ben(m1,m2 -> 0) == massless.
   Sends BOTH masses to zero; every 2-mass PaVe -> its massless limit
   (Bf/Cf1/Cf2/boxF), and the assembled amplitude must reproduce
   amplitudeMass0Stefan (notebook-validated) order by order in eps.
   Validates the freshly-generated 2-mass reduction ampRaw_2mass.txt end-to-end.
   Run:  wolframscript -file capstone_m0_2mass.wl
============================================================================ *)
ep=Global`Epsilon; mu2=1; ORD=3;
cG=cGsym; cGRule=cGsym->1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
MSBarFac=(16 Pi^2/I)(Exp[EulerGamma]/(4 Pi))^ep;
hypE[zz_]:=Normal[N[Series[Hypergeometric2F1[-ep,-ep,1-ep,N[zz,30]],{ep,0,ORD+3}],25]];
boxF[a_,b_,q_]:=Module[{cc=q-a-b}, 2 I cG/ep^2 (1/(a b))(
   (-(b+cc)/(a (b/mu2)))^ep hypE[1-b/(b+cc)]
 + (-(a+cc)/(a (b/mu2)))^ep hypE[1-a/(a+cc)]
 - (-(a+cc)(b+cc)/((a+b+cc) a (b/mu2)))^ep hypE[1-(a b)/((a+cc)(b+cc))])];
Bf[x_]:=I cG/(ep(1-2 ep))(-mu2/x)^ep;
Cf1[x_]:=I cG/ep^2 (1/x)(-mu2/x)^ep;
Cf2[x_,y_]:=I cG/ep^2 (1/(x-y))((-mu2/x)^ep-(-mu2/y)^ep);

mapML2[sv_,tv_,uv_]:=With[{Q2=sv+tv+uv},{
  B0[m1^2,0,m1^2]->0, B0[m2^2,0,m2^2]->0,
  B0[s,m1^2,m2^2]:>Bf[sv], B0[t,0,m1^2]:>Bf[tv], B0[u,0,m2^2]:>Bf[uv],
  B0[-m1^2-m2^2+s+t+u,m1^2,m2^2]:>Bf[Q2],
  C0[0,m1^2,t,m1^2,m1^2,0]:>Cf1[tv], C0[0,m2^2,u,m2^2,m2^2,0]:>Cf1[uv],
  C0[0,s,-m1^2-m2^2+s+t+u,m1^2,m1^2,m2^2]:>Cf2[sv,Q2],
  C0[0,s,-m1^2-m2^2+s+t+u,m2^2,m2^2,m1^2]:>Cf2[sv,Q2],
  C0[m1^2,m2^2,s,m1^2,0,m2^2]:>Cf1[sv],
  C0[m1^2,u,-m1^2-m2^2+s+t+u,m1^2,0,m2^2]:>Cf2[uv,Q2],
  C0[m2^2,t,-m1^2-m2^2+s+t+u,m2^2,0,m1^2]:>Cf2[tv,Q2],
  D0[0,u,m1^2,s,m2^2,-m1^2-m2^2+s+t+u,m2^2,m2^2,0,m1^2]:>boxF[uv,sv,Q2],
  D0[m1^2,t,-m1^2-m2^2+s+t+u,s,0,m2^2,m1^2,0,m1^2,m2^2]:>boxF[sv,tv,Q2]}];

$dir=DirectoryName[$InputFileName]; $f2=ParentDirectory[$dir];
$f0=FileNameJoin[{ParentDirectory[$f2],"for_stefan"}];
raw=Import[FileNameJoin[{$f2,"ampRaw_2mass.txt"}],"Text"];
amp0raw=ToExpression[raw,InputForm]/.{SUNN->3,SDF[__]->1,SUNTF[__]->1,SMP["g_s"]->1,
   ca->1,ci->1,cj->1,ck->1,cl->1};
ev[sv_,tv_,uv_]:=Module[{a},
  a=amp0raw/.mapML2[sv,tv,uv]/.{m1->0,m2->0}/.{s->sv,t->tv,u->uv};
  a=a MSBarFac/.cGRule;
  Table[N[SeriesCoefficient[Series[a,{ep,0,ORD}],k],10],{k,-2,ORD}]];

Get[FileNameJoin[{$f0,"amplitude_functions.wl"}]];
Print["pt {s,t,u}:  ratio amp2(m1,m2->0)/amp0Stefan  (eps^-2..eps^3)"];
Do[Module[{sv=p[[1]],tv=p[[2]],uv=p[[3]],a2,a0},
  a2=ev[sv,tv,uv]; a0=amplitudeMass0Stefan[sv,tv,uv,sv+tv+uv];
  Print[p,"  ",Table[If[Abs[a0[[k]]]>10^-12,N[Re[a2[[k]]/a0[[k]]],8],"-"],{k,1,6}]]],
 {p,{{-1,-1,-1},{-2,-1,-1},{-1,-2,-1/2},{-3,-1,-1/2},{-1,-5,-2}}}];
