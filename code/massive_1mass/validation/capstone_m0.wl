(* ============================================================================
   CAPSTONE VALIDATION:  amplitudeMass1Ben(m^2 -> 0)  ==  amplitudeMass0Stefan
   ----------------------------------------------------------------------------
   Validates the freshly-generated FeynCalc MASSIVE reduction (ampRaw_1mass.txt)
   against the already-Stefan-validated massless amplitude, by taking m^2 -> 0:
   every massive PaVe collapses to its massless limit (BFunc/CFunc/DFunc forms),
   and the assembled amplitude must reproduce amplitudeMass0Stefan order by order.

   RESULT:  ratio amp1(m->0)/amp0 = 1.000000 at every order eps^-2..eps^3, at
   5 kinematic points (2 symmetric + 3 fully asymmetric).  => the massive
   reduction + the master->massless-limit map (incl. the box arg assignment) are
   correct.

   Box arg map (pinned by this test):
     D0_a = D0(0,u,mq2,s; mq2,Q2; mq2,mq2,0,mq2)  ->  boxF[u, s, Q2]
     D0_b = D0(mq2,t,Q2,s; 0,mq2; mq2,0,mq2,mq2)  ->  boxF[s, t, Q2]

   Run:  wolframscript -file capstone_m0.wl
============================================================================ *)
ep=Global`Epsilon; mu2=1; ORD=3;
cG=cGsym;
cGRule=cGsym->1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
MSBarFac=(16 Pi^2/I)(Exp[EulerGamma]/(4 Pi))^ep;

(* massless master closed forms (cG / i / mu2 convention; == package boxF/StefanRules) *)
hypE[zz_]:=Normal[N[Series[Hypergeometric2F1[-ep,-ep,1-ep,N[zz,30]],{ep,0,ORD+3}],25]];
boxF[a_,b_,q_]:=Module[{cc=q-a-b}, 2 I cG/ep^2 (1/(a b))(
   (-(b+cc)/(a (b/mu2)))^ep hypE[1-b/(b+cc)]
 + (-(a+cc)/(a (b/mu2)))^ep hypE[1-a/(a+cc)]
 - (-(a+cc)(b+cc)/((a+b+cc) a (b/mu2)))^ep hypE[1-(a b)/((a+cc)(b+cc))])];
Bf[x_]:=I cG/(ep(1-2 ep))(-mu2/x)^ep;                               (* B0(x;0,0) *)
Cf1[x_]:=I cG/ep^2 (1/x)(-mu2/x)^ep;                                (* C0(0,0,x;0,0,0) *)
Cf2[x_,y_]:=I cG/ep^2 (1/(x-y))((-mu2/x)^ep-(-mu2/y)^ep);           (* C0(0,x,y;0,0,0) *)

(* massive PaVe -> massless limit (mq^2 -> 0). boxA,boxB are the box values. *)
mapML[bA_,bB_]:={B0[mq^2,0,mq^2]->0,
  B0[s,mq^2,mq^2]:>Bf[ss],B0[t,0,mq^2]:>Bf[tt],B0[u,0,mq^2]:>Bf[uu],
  B0[-2 mq^2+s+t+u,mq^2,mq^2]:>Bf[ss+tt+uu],
  C0[0,mq^2,t,mq^2,mq^2,0]:>Cf1[tt],C0[0,mq^2,u,mq^2,mq^2,0]:>Cf1[uu],
  C0[mq^2,mq^2,s,mq^2,0,mq^2]:>Cf1[ss],
  C0[0,s,-2 mq^2+s+t+u,mq^2,mq^2,mq^2]:>Cf2[ss,ss+tt+uu],
  C0[mq^2,t,-2 mq^2+s+t+u,mq^2,0,mq^2]:>Cf2[tt,ss+tt+uu],
  C0[mq^2,u,-2 mq^2+s+t+u,mq^2,0,mq^2]:>Cf2[uu,ss+tt+uu],
  D0[0,u,mq^2,s,mq^2,-2 mq^2+s+t+u,mq^2,mq^2,0,mq^2]:>bA,
  D0[mq^2,t,-2 mq^2+s+t+u,s,0,mq^2,mq^2,0,mq^2,mq^2]:>bB};

$dir=DirectoryName[$InputFileName]; $f1=ParentDirectory[$dir];
$f0=FileNameJoin[{ParentDirectory[$f1],"for_stefan"}];
raw=Import[FileNameJoin[{$f1,"ampRaw_1mass.txt"}],"Text"];
amp0raw=ToExpression[raw,InputForm]/.{SUNN->3,SDF[__]->1,SUNTF[__]->1,SMP["g_s"]->1,
   ca->1,ci->1,cj->1,ck->1,cl->1};
ev[sv_,tv_,uv_]:=Module[{a,qv=sv+tv+uv,A,B},
  A=boxF[uv,sv,qv]; B=boxF[sv,tv,qv];
  a=amp0raw/.mapML[A,B]/.mq->0/.{s->sv,t->tv,u->uv}/.{ss->sv,tt->tv,uu->uv};
  a=a MSBarFac/.cGRule;
  Table[N[SeriesCoefficient[Series[a,{ep,0,ORD}],k],10],{k,-2,ORD}]];

Get[FileNameJoin[{$f0,"amplitude_functions.wl"}]];   (* amplitudeMass0Stefan reference *)
Print["pt {s,t,u} :  ratio amp1(m->0)/amp0Stefan  (eps^-2 .. eps^3)"];
Do[Module[{sv=p[[1]],tv=p[[2]],uv=p[[3]],a1,a0},
  a1=ev[sv,tv,uv]; a0=amplitudeMass0Stefan[sv,tv,uv,sv+tv+uv];
  Print[p,"  ",Table[If[Abs[a0[[k]]]>10^-12,N[Re[a1[[k]]/a0[[k]]],8],"-"],{k,1,6}]]],
 {p,{{-1,-1,-1},{-2,-1,-1},{-1,-2,-1/2},{-3,-1,-1/2},{-1,-5,-2}}}];
