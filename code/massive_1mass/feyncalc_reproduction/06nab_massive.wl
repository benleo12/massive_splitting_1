(* 06nab: MASSIVE non-abelian radiator — the correct massive analogue of
   for_stefan/feyncalc_reproduction/01+02 (NOT the diagsl2 full-colour run).
   Process : {S[1,{ck,cj}]} -> {S[13,{1,ci}], V[50,{ca}], -S[13,{1,cl}]}   (anti-squark!)
   Diagrams: 14, (13+15), (5+12)  with the (I/2) SUNF colour projections
   Massive : SQCDBGF_massive model (mq), external squark shells k1^2=k3^2=mq^2
   NOTE: no TrickMandelstam (5 scales {s,t,u,Q^2,mq^2} break the Mandelstam relation,
         since s+t+u = Q^2 + 2 mq^2); keep (s,t,u,mq^2) independent.
   Output : ampRaw_nab_1mass.txt (k1.eps coefficient, colour-stripped by nabprefac)
            + the master list.  m->0 must reproduce Hoeche's MyBox[s,t,u]. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];

diags2=InsertFields[CreateTopologies[1,1->3,ExcludeTopologies->{Tadpoles}],
  {S[1,{ck,cj}]}->{S[13,{1,ci}],V[50,{ca}],-S[13,{1,cl}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],F[_]}];
conv[idx_,lm_]:=FCFAConvert[CreateFeynAmp[DiagramExtract[diags2,{idx}]],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->False,LoopMomenta->lm,SMP->True,Contract->True,
  DropSumOver->True];   (* NO MQU->0 : keep the squark mass *)
polk2=Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];
tid[a_]:=TID[FCE[a],l,UsePaVeBasis->True,ToPaVe->True]/.{polk2->0};

amp12$0 = conv[14,{l}];
amps1113$0 = Simplify[ScalarProductExpand[
  (conv[13,{-l+k1+k3}] /. {SUNTF[{a_},i_,j_] SUNTF[{a_},m_,l_] SUNTF[{b_},k_,m_] :>
       (I/2) SUNF[b,a,SUNIndex[Glu6]] SUNTF[{SUNIndex[Glu6]},k,l] SUNTF[{a},i,j]})
 +(conv[15,{l+k2}] /. {SUNTF[{a_},i_,j_] SUNTF[{b_},m_,l_] SUNTF[{a_},k_,m_] :>
       (I/2) SUNF[SUNIndex[Glu6],b,a] SUNTF[{a},k,l] SUNTF[{SUNIndex[Glu6]},i,j]})
  ] /. {polk2->0}];
amps510$0 = {
  SUNSimplify[conv[5,{l}]] /. {SUNTF[{a_,b_},i_,j_] :>
       (I/2) SUNF[a,b,SUNIndex[Glu6]] SUNTF[{SUNIndex[Glu6]},i,j]} /. {Glu5->Glu6,Glu6->Glu5},
  SUNSimplify[conv[12,{l}]] /. {SUNTF[{a_,b_},i_,j_] :>
       (I/2) SUNF[a,b,SUNIndex[Glu6]] SUNTF[{SUNIndex[Glu6]},i,j]} };
Print["LCs: ", LeafCount/@{amp12$0, amps1113$0, Total[amps510$0]}];

Print["TID..."];
amp12$1 = tid[amp12$0]; amps1113$1 = tid[amps1113$0]; amps510$1 = tid[Total[amps510$0]];
Print["TID done."];

nabprefac = (1/(16 Pi^2)) (SMP["g_s"]^3 SUNF[SUNIndex[ca],SUNIndex[Glu5],SUNIndex[Glu6]]*
   SUNTF[{SUNIndex[Glu5]},SUNFIndex[ci],SUNFIndex[cj]] SUNTF[{SUNIndex[Glu6]},SUNFIndex[ck],SUNFIndex[cl]]);
red[a_] := PaVeReduce[FeynAmpDenominatorExplicit[a], PaVeAutoReduce->True];
Print["PaVeReduce 1/3..."]; a12   = red[amp12$1];   Print["  LC ",LeafCount[a12]];
Print["PaVeReduce 2/3..."]; a1113 = red[amps1113$1];Print["  LC ",LeafCount[a1113]];
Print["PaVeReduce 3/3..."]; a510  = red[amps510$1]; Print["  LC ",LeafCount[a510]];

sum = a12 + a1113 + a510;
nabRaw = SUNSimplify[sum/nabprefac, SUNNToCACF->False] // Simplify;
polk1=Pair[Momentum[k1,D],Momentum[Polarization[k2,-I],D]];
polk3=Pair[Momentum[k3,D],Momentum[Polarization[k2,-I],D]];
coeffK1 = Coefficient[nabRaw /. {polk3->0}, polk1];
Print["MASSIVE non-abelian k1-coeff LC: ", LeafCount[coeffK1]];
Print["mq present? ", !FreeQ[coeffK1, mq]];
Print["A0 args: ", Union[Cases[coeffK1,A0[x__]:>{x},Infinity]]];
Print["B0 args: ", Union[Cases[coeffK1,B0[x__]:>{x},Infinity]]];
Print["C0 args: ", Union[Cases[coeffK1,C0[x__]:>{x},Infinity]]];
Print["D0 args: ", Union[Cases[coeffK1,D0[x__]:>{x},Infinity]]];
Export[FileNameJoin[{base,"for_stefan_1mass","ampRaw_nab_1mass.txt"}],
   ToString[FCE[coeffK1], InputForm], "Text"];
Print["saved ampRaw_nab_1mass.txt"];
(* m->0 sanity: master call set should collapse to the massless non-abelian set *)
Print["--- m->0 master call set ---"];
c0 = coeffK1 /. mq->0;
Print["B0: ", Union[Cases[c0,B0[x__]:>{x},Infinity]]];
Print["C0: ", Union[Cases[c0,C0[x__]:>{x},Infinity]]];
Print["D0: ", Union[Cases[c0,D0[x__]:>{x},Infinity]]];
