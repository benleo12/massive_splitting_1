PrependTo[$Path,"/tmp/fc93"]; $LoadAddOns={"FeynArts","FeynHelpers"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0;
SetDirectory["/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop/collaborator_massless"];
Get["/tmp/ab_coeffK1.mx"];   (* restores abK1 *)
ep=Epsilon; mu2=1; cG=Subscript[c,"\[CapitalGamma]"];
cGRule = cG -> 1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
hyp[z_]:=1+ep^2 PolyLog[2,z];
Bub[x_]:=I cG/(ep(1-2 ep))(-mu2/x)^ep;
Tri[x_]:=I cG/ep^2 (1/x)(-mu2/x)^ep;
Tri2[x_,y_]:=I cG/ep^2 (1/(x-y))((-mu2/x)^ep-(-mu2/y)^ep);
BoxF[a_,b_,d_]:=2 I cG/ep^2 (1/(a b))((-(b+d)/(a(b/mu2)))^ep hyp[1-b/(b+d)]
 +(-(a+d)/(a(b/mu2)))^ep hyp[1-a/(a+d)]
 -(-(a+d)(b+d)/((a+b+d)a(b/mu2)))^ep hyp[1-(a b)/((a+d)(b+d))]);

(* Cell 90 MyBoxAmpAbLS1, colour structures C1,C2; colour-STRIPPED both ->2/3 *)
cstrip = 2/3;  (* (1/2)(1+1) - 1/SUNN with SUNN->3 *)
MyAbLS1[s_,t_,u_] := cstrip (
   (I/t)((-u)((Bub[s+t+u]-Bub[s])/(s+t+u-s)) - 2(s+t+u-t)Tri2[s+t+u,t] - Bub[s] - s t BoxF[s,t,u])
 + (I/t)( t ((Bub[s+t+u]-Bub[s])/(s+t+u-s)) + 2(s+t+u-u)Tri2[s+t+u,u] - 2 s Tri[s] + u s BoxF[u,s,t]) );
(* full reference: MyBoxAmpAbLS1 * (I/(16 Pi^2)) * g_s^3, g_s->1 *)
stefanRef[s_,t_,u_] := MyAbLS1[s,t,u] (I/(16 Pi^2));

(* our abK1: strip colour + coupling, evaluate masters via PaX *)
pfx=(Exp[-EulerGamma]/(4 Pi))^(-Epsilon);
ourStrip = abK1 /. {SUNN->3, SDF[__]->1, SUNFDelta[__]->1, SUNTF[__]->1, SMP["g_s"]->1};
paxVal[e_]:=e/.{B0[x__]:>PaXEvaluate[B0[x],PaXImplicitPrefactor->pfx],
   C0[x__]:>PaXEvaluate[C0[x],PaXImplicitPrefactor->pfx],
   D0[x__]:>PaXEvaluate[D0[x],PaXImplicitPrefactor->pfx]};
ser[e_]:=Series[e/.ScaleMu->1,{Epsilon,0,0}];
ourS[pt_]:=Module[{sv,tv,uv,qv},{sv,tv,uv}=pt;qv=sv+tv+uv;
   ser[paxVal[ourStrip/.{s->sv,t->tv,u->uv,Q^2->qv,Q->Sqrt[qv]}]]];
stfS[pt_]:=Module[{sv,tv,uv},{sv,tv,uv}=pt; ser[stefanRef[sv,tv,uv]/.cGRule]];

pts={{-1,-1,-1},{-2,-1,-1},{-1,-2,-1/2},{-3,-2,-2}};
N1=Series[ourS[pts[[1]]]/stfS[pts[[1]]],{Epsilon,0,2}];
Print["N(eps) from P1: ", N[Normal[N1],8], "   (leading -16Pi^2 I = ",N[-16 Pi^2 I,8],")"];
Do[Module[{resid},
  resid = Normal[ourS[pt]-N1 stfS[pt]];
  resid = Table[N[SeriesCoefficient[Series[resid,{Epsilon,0,0}],k],8],{k,-2,0}];
  Print["pt ",pt,":  our - N*Stefan = ",resid,"   max|.| ",Max[Abs[resid]]];
],{pt,pts}];
