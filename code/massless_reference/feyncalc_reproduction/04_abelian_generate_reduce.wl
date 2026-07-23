PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
collab="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop/collaborator_massless";
SetDirectory[collab]; modelAbs=FileNameJoin[{collab,"SQCDBGF","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,0,0];
CF=(SUNN^2-1)/(2 SUNN); CA=SUNN;
(* abelian color rules *)
TaTAbaRule := {SUNTF[{a_},i_,j_](SUNTF[{b_,a_},k_,l_]+SUNTF[{a_,b_},k_,l_]) ->
   SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2)+(-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]};
proj1113 := {SUNTF[{a_},k_,m_]SUNTF[{a_},i_,j_]SUNTF[{b_},m_,l_] ->
   (1/2)((-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]+SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2))};
proj510 := {SUNTF[{a_},i_,j_]SUNTF[{b_,a_},k_,l_] ->
   (1/2)((-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]+SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2))};

diagsl2=InsertFields[CreateTopologies[1,1->3,ExcludeTopologies->{Tadpoles}],
  {S[1,{cl,cj}]}->{S[13,{1,ci}],V[50,{ca}],S[13,{1,ck}]}, InsertionLevel->{Classes},
  Model->modelAbs, GenericModel->modelAbs, ExcludeParticles->{S[1],S[14],F[_]}];
conv[idx_,lm_]:=FCFAConvert[CreateFeynAmp[DiagramExtract[diagsl2,{idx}]],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->False,LoopMomenta->lm,SMP->True,Contract->True,
  DropSumOver->True,FinalSubstitutions->{MQU[__]->0}];
polk2=Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];

(* Cell 81 *)
ampsl12$0 = (conv[1,{l}] /. TaTAbaRule) + ((conv[2,{l}] /. {cj->cl,cl->cj}) /. TaTAbaRule);
(* Cell 79 *)
ampsal1113$0 = Simplify[ScalarProductExpand[{
   conv[9,{-l+k1+k3}] /. proj1113, conv[11,{l+k2}] /. proj1113 }] /. {polk2->0}];
(* Cell 84 *)
ampsal510$0 = { SUNSimplify[conv[5,{l}]] /. {cl->cj,cj->cl} /. proj510,
                SUNSimplify[conv[8,{l}]] /. proj510 };
Print["amp[0] LCs: ", LeafCount/@{ampsl12$0, Total@Flatten@{ampsal1113$0}, Total@Flatten@{ampsal510$0}}];

tid[a_]:=TID[FCE[a],l,UsePaVeBasis->True,ToPaVe->True]/.{polk2->0};
red[a_]:=Module[{r}, r=a/.(h:A0|B0|C0|D0|PaVe)[x__]:>TrickMandelstam[h[x],{s,t,u,Q^2}];
   r=FeynAmpDenominatorExplicit[r]; r=PaVeReduce[r,PaVeAutoReduce->True];
   r/.(h:A0|B0|C0|D0|PaVe)[x__]:>TrickMandelstam[h[x],{s,t,u,Q^2}]];

Print["TID+reduce ampsl12..."];   a12 = red[tid[ampsl12$0]];
Print["TID+reduce ampsal1113..."];a1113 = red[tid[Total@Flatten@{ampsal1113$0}]];
Print["TID+reduce ampsal510..."]; a510 = red[tid[Total@Flatten@{ampsal510$0}]];

(* Cell 87 assembly: k1 projection = set Pair[k3,pol]->0, coeff of Pair[k1,pol] *)
sum = a12 + a1113 + a510;
polk1=Pair[Momentum[k1,D],Momentum[Polarization[k2,-I],D]];
polk3=Pair[Momentum[k3,D],Momentum[Polarization[k2,-I],D]];
abK1 = Coefficient[sum /. {polk3->0}, polk1];
Print["abelian k1-coeff LC: ", LeafCount[abK1]];
Print["color heads: ", Union[Head/@Cases[abK1,(SUNTF|SUNFDelta|SUNF)[___],Infinity]]];
DumpSave["/tmp/ab_coeffK1.mx", abK1];
Print["saved /tmp/ab_coeffK1.mx"];

(* diagnostics *)
Print["--- diag ---"];
Print["Polarization present in sum? ", !FreeQ[sum, Polarization]];
Print["distinct Pair[...Polarization...] in sum: ",
  Take[#, Min[6,Length[#]]]&@Union[Cases[sum, Pair[a_,b_]/;!FreeQ[{a,b},Polarization]:>Pair[a,b], Infinity]]];
DumpSave["/tmp/ab_sum.mx", sum];
