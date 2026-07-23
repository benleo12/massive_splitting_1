(* 05n: full-colour reduction of the complete 21-diagram physical amplitude.
   Loads the converted amps from /tmp/full_amps_1mass.mx, keeps FULL colour,
   projects onto k1.eps(k2) (and separately k3), TID + PaVe reduce -> complete
   Fig-14 radiator (abelian + non-abelian).  Masters == the validated abelian set. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base]; SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
Get["/tmp/full_amps_1mass.mx"];   (* -> amps (list of 21) *)
Print["loaded ", Length[amps], " amps"];

polk2=Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];
tid[a_]:=TID[FCE[a],l,UsePaVeBasis->True,ToPaVe->True]/.{polk2->0};
red[a_]:=PaVeReduce[FeynAmpDenominatorExplicit[a],PaVeAutoReduce->True];

(* reduce diagram-by-diagram (robust) then sum *)
Print["reducing ", Length[amps], " diagrams..."];
redAmps = Table[
   Module[{r=Quiet[red[tid[amps[[i]]]]]}, Print["  diag ",i," done, LC=",LeafCount[r]]; r],
   {i, Length[amps]}];
sum = FCE[Total[redAmps]];
polk1=SPD[k1,Polarization[k2,-I]]; polk3=SPD[k3,Polarization[k2,-I]];
full = sum;
abK1 = Coefficient[sum /. {polk3->0}, polk1];
abK3 = Coefficient[sum /. {polk1->0}, polk3];
Print["full k1-coeff LC: ", LeafCount[abK1], "   non-abelian present? ", !FreeQ[abK1,SUNF]||!FreeQ[abK1,SUNN^2]];
Print["colour heads: ", Union[Head/@Cases[abK1,(SUNTF|SUNFDelta|SUNF)[___],Infinity]]];
out=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","ampRaw_full_1mass.txt"}];
Export[out, ToString[abK1,InputForm], "Text"];
Export[StringReplace[out,"full_1mass"->"full_k3_1mass"], ToString[abK3,InputForm], "Text"];
Print["saved ", out];
DumpSave["/tmp/ampRaw_full_1mass.mx", {abK1,abK3}];
