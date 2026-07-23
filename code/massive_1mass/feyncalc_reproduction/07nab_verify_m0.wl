(* 07nab: VALIDATE the massive non-abelian radiator by its m->0 limit against
   Hoeche's hand-derived closed form MyBox[s,t,u] (the massless non-abelian).
   Mirrors for_stefan/feyncalc_reproduction/03_nonabelian_verify.wl, but the input is
   OUR massive amplitude (ampRaw_nab_1mass.txt) with mq -> 0.
   Success = the ratio our/MyBox is a kinematics-independent constant N(eps)
   (leading -16 Pi^2 I), identical at every point. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts","FeynHelpers"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[FileNameJoin[{base,"collaborator_massless"}]];

coeffK1m = ToExpression[Import[FileNameJoin[{base,"for_stefan_1mass","ampRaw_nab_1mass.txt"}],"Text"], InputForm];
coeffK1  = FCI[coeffK1m /. mq -> 0];        (* massless limit of OUR massive result *)
Print["our m->0 amplitude LC: ", LeafCount[coeffK1]];

ep=Epsilon; mu2=1; cG=Subscript[c,"\[CapitalGamma]"];
cGRule = cG -> 1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
hyp[z_]:=1+ep^2 PolyLog[2,z];
Bub[x_]:=I cG/(ep(1-2 ep))(-mu2/x)^ep;
Tri[x_]:=I cG/ep^2 (1/x)(-mu2/x)^ep;
BoxF[a_,b_,d_]:=2 I cG/ep^2 (1/(a b))((-(b+d)/(a(b/mu2)))^ep hyp[1-b/(b+d)]
 +(-(a+d)/(a(b/mu2)))^ep hyp[1-a/(a+d)]
 -(-(a+d)(b+d)/((a+b+d)a(b/mu2)))^ep hyp[1-(a b)/((a+d)(b+d))]);
(* Hoeche's hand-derived massless non-abelian reference (cell 40) *)
MyBox[s_,t_,u_]:=(I/t)(Bub[s+t+u]+2 t Tri[t]+2 u Tri[u]+2 s Tri[s]
   -s t BoxF[s,t,u]-u s BoxF[u,s,t]+2 t u BoxF[t,u,s]);

pfx=(Exp[-EulerGamma]/(4 Pi))^(-Epsilon);
paxVal[e_]:=e/.{B0[x__]:>PaXEvaluate[B0[x],PaXImplicitPrefactor->pfx],
   C0[x__]:>PaXEvaluate[C0[x],PaXImplicitPrefactor->pfx],
   D0[x__]:>PaXEvaluate[D0[x],PaXImplicitPrefactor->pfx]};
ser[e_]:=Series[e/.ScaleMu->1,{Epsilon,0,0}];
ourS[pt_]:=Module[{sv,tv,uv,qv},{sv,tv,uv}=pt; qv=sv+tv+uv;
   ser[paxVal[coeffK1/.{s->sv,t->tv,u->uv,Q^2->qv,Q->Sqrt[qv]}]]];
stfS[pt_]:=Module[{sv,tv,uv},{sv,tv,uv}=pt; ser[MyBox[sv,tv,uv]/.cGRule]];

pts={{-1,-1,-1},{-2,-1,-1},{-1,-2,-1/2},{-3,-2,-2}};
Do[Module[{r},
   r = Series[ourS[p]/stfS[p], {Epsilon,0,2}];
   Print["pt ",p,"  our(m->0)/MyBox = ", N[Normal[r],8]]],
 {p,pts}];
Print["  (should be the SAME constant at every point; leading = -16 Pi^2 I = ", N[-16 Pi^2 I,8],")"];
