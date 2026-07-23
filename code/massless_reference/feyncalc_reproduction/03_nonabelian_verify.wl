PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts","FeynHelpers"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0;
SetDirectory["/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop/collaborator_massless"];
Get["/tmp/nab_coeffK1.mx"];
ep=Epsilon; mu2=1; cG=Subscript[c,"\[CapitalGamma]"];
cGRule = cG -> 1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
hyp[z_]:=1+ep^2 PolyLog[2,z];
Bub[x_]:=I cG/(ep(1-2 ep))(-mu2/x)^ep;
Tri[x_]:=I cG/ep^2 (1/x)(-mu2/x)^ep;
BoxF[a_,b_,d_]:=2 I cG/ep^2 (1/(a b))((-(b+d)/(a(b/mu2)))^ep hyp[1-b/(b+d)]
 +(-(a+d)/(a(b/mu2)))^ep hyp[1-a/(a+d)]
 -(-(a+d)(b+d)/((a+b+d)a(b/mu2)))^ep hyp[1-(a b)/((a+d)(b+d))]);
MyBox[s_,t_,u_]:=(I/t)(Bub[s+t+u]+2 t Tri[t]+2 u Tri[u]+2 s Tri[s]
   -s t BoxF[s,t,u]-u s BoxF[u,s,t]+2 t u BoxF[t,u,s]);
pfx=(Exp[-EulerGamma]/(4 Pi))^(-Epsilon);
paxVal[e_]:=e/.{B0[x__]:>PaXEvaluate[B0[x],PaXImplicitPrefactor->pfx],
   C0[x__]:>PaXEvaluate[C0[x],PaXImplicitPrefactor->pfx],
   D0[x__]:>PaXEvaluate[D0[x],PaXImplicitPrefactor->pfx]};
(* return Laurent series object eps^-2..eps^0 *)
ser[e_]:=Series[e/.ScaleMu->1,{Epsilon,0,0}];

ourS[pt_]:=Module[{sv,tv,uv,qv},{sv,tv,uv}=pt; qv=sv+tv+uv;
   ser[paxVal[coeffK1/.{s->sv,t->tv,u->uv,Q^2->qv,Q->Sqrt[qv]}]]];
stfS[pt_]:=Module[{sv,tv,uv},{sv,tv,uv}=pt; ser[MyBox[sv,tv,uv]/.cGRule]];

pts={{-1,-1,-1},{-2,-1,-1},{-1,-2,-1/2},{-3,-2,-2}};
(* N(eps) from P1 as series ratio *)
N1 = Series[ourS[pts[[1]]]/stfS[pts[[1]]], {Epsilon,0,2}];
Print["N(eps) from P1 (PaX/Stefan normalisation): ", N[Normal[N1],8]];
Print["  (leading should be -16 Pi^2 I = ", N[-16 Pi^2 I,8], ")"];
Print[""];
Do[Module[{resid},
   resid = Normal[ourS[pt] - N1 stfS[pt]];
   resid = Table[N[SeriesCoefficient[Series[resid,{Epsilon,0,0}],k],8],{k,-2,0}];
   Print["pt ",pt,":  our - N*Stefan  per eps = ", resid, "   max|.| ", Max[Abs[resid]]];
],{pt,pts}];
