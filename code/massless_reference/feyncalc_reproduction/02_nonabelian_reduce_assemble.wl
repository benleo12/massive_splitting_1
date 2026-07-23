PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
collab="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop/collaborator_massless";
SetDirectory[collab];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,0,0];
{amp12$1, amps1113$1, amps510$1} = Get["/tmp/nab_amp1.mx"];

nabprefac = (1/(16 Pi^2)) (SMP["g_s"]^3 SUNF[SUNIndex[ca],SUNIndex[Glu5],SUNIndex[Glu6]]*
   SUNTF[{SUNIndex[Glu5]},SUNFIndex[ci],SUNFIndex[cj]] SUNTF[{SUNIndex[Glu6]},SUNFIndex[ck],SUNFIndex[cl]]);

(* amp[2]: reduce PaVe tensors -> scalars via FeynCalc PaVeReduce, TrickMandelstam args *)
red[a_] := Module[{r},
   r = a /. (h:A0|B0|C0|D0|PaVe)[x__] :> TrickMandelstam[h[x],{s,t,u,Q^2}];
   r = FeynAmpDenominatorExplicit[r];
   r = PaVeReduce[r, PaVeAutoReduce->True];
   r /. (h:A0|B0|C0|D0|PaVe)[x__] :> TrickMandelstam[h[x],{s,t,u,Q^2}] ];

Print["PaVeReduce amp12..."];  a12 = red[amp12$1];
Print["  done LC ",LeafCount[a12]];
Print["PaVeReduce amps1113..."]; a1113 = red[amps1113$1];
Print["  done LC ",LeafCount[a1113]];
Print["PaVeReduce amps510..."]; a510 = red[amps510$1];
Print["  done LC ",LeafCount[a510]];

sum = a12 + a1113 + a510;
(* divide by nabprefac, reduce color to scalar *)
nabRaw = SUNSimplify[sum/nabprefac, SUNNToCACF->False] // Simplify;
(* cell 36: set Pair[k3,pol]->0, keep Pair[k1,pol] coefficient *)
polk = Pair[Momentum[k1,D],Momentum[Polarization[k2,-I],D]];
nab1 = nabRaw /. {Pair[Momentum[k3,D],Momentum[Polarization[k2,-I],D]]->0};
coeffK1 = Coefficient[nab1, polk];
Print["coeff of k1.eps : LC ", LeafCount[coeffK1]];
Print["distinct scalar masters in our reduced amp:"];
Print["  B0 args: ", Union[Cases[coeffK1, B0[x__]:>{x}, Infinity]]];
Print["  C0 args: ", Union[Cases[coeffK1, C0[x__]:>{x}, Infinity]]];
Print["  D0 args: ", Union[Cases[coeffK1, D0[x__]:>{x}, Infinity]]];
DumpSave["/tmp/nab_coeffK1.mx", coeffK1];
Print["saved /tmp/nab_coeffK1.mx"];
