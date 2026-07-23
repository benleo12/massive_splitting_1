(* 04m: massive-squark version of for_stefan/feyncalc_reproduction/04_abelian_generate_reduce.wl
   SAME 6 abelian diagrams {1,2,9,11,5,8} and SAME colour projections, but:
     - load the local SQCDBGF_massive model (gen-1 up-squark mass kept symbolic = mq),
     - NO  MQU[__]->0  substitution,
     - external squark legs on massive shell: SetMandelstam[...,Q,0,mq,mq]
       (k2 = gluon massless; k1,k3 = squarks, mass mq) -> s+t+u = Q^2 + 2 mq^2.
   Output: abelian k1-coeff (massive) saved to /tmp/ab_coeffK1_massive.mx, master list,
   and a m->0 sanity check of the master CALL SET vs the massless masters. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
CF=(SUNN^2-1)/(2 SUNN); CA=SUNN;
TaTAbaRule := {SUNTF[{a_},i_,j_](SUNTF[{b_,a_},k_,l_]+SUNTF[{a_,b_},k_,l_]) ->
   SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2)+(-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]};
proj1113 := {SUNTF[{a_},k_,m_]SUNTF[{a_},i_,j_]SUNTF[{b_},m_,l_] ->
   (1/2)((-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]+SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2))};
proj510 := {SUNTF[{a_},i_,j_]SUNTF[{b_,a_},k_,l_] ->
   (1/2)((-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]+SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2))};

diagsl2=InsertFields[CreateTopologies[1,1->3,ExcludeTopologies->{Tadpoles}],
  {S[1,{cl,cj}]}->{S[13,{1,ci}],V[50,{ca}],S[13,{1,ck}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],F[_]}];
(* NOTE: NO MQU->0 in FinalSubstitutions -> the squark mass mq is kept *)
conv[idx_,lm_]:=FCFAConvert[CreateFeynAmp[DiagramExtract[diagsl2,{idx}]],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->False,LoopMomenta->lm,SMP->True,Contract->True,DropSumOver->True];
polk2=Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];

ampsl12$0   = (conv[1,{l}] /. TaTAbaRule) + ((conv[2,{l}] /. {cj->cl,cl->cj}) /. TaTAbaRule);
ampsal1113$0= Simplify[ScalarProductExpand[{
   conv[9,{-l+k1+k3}] /. proj1113, conv[11,{l+k2}] /. proj1113 }] /. {polk2->0}];
ampsal510$0 = { SUNSimplify[conv[5,{l}]] /. {cl->cj,cj->cl} /. proj510,
                SUNSimplify[conv[8,{l}]] /. proj510 };
Print["amp[0] LCs: ", LeafCount/@{ampsl12$0, Total@Flatten@{ampsal1113$0}, Total@Flatten@{ampsal510$0}}];

tid[a_]:=TID[FCE[a],l,UsePaVeBasis->True,ToPaVe->True]/.{polk2->0};
(* NO TrickMandelstam: the 5-scale list {s,t,u,Q^2,mq^2} breaks the Mandelstam
   relation (now s+t+u = Q^2 + 2 mq^2) and zeroes the result.  Keep (s,t,u,mq^2)
   independent; Q^2 = s+t+u-2 mq^2 is reconstructed downstream. *)
red[a_]:=PaVeReduce[FeynAmpDenominatorExplicit[a],PaVeAutoReduce->True];

Print["TID+reduce ampsl12..."];   a12  = red[tid[ampsl12$0]];
Print["TID+reduce ampsal1113..."];a1113= red[tid[Total@Flatten@{ampsal1113$0}]];
Print["TID+reduce ampsal510..."]; a510 = red[tid[Total@Flatten@{ampsal510$0}]];

(* extract in FCE/SPD form (PaVeReduce output is FCI; Pair[Momentum[k1,D],..]
   no longer matches -- use SPD[k1,Polarization[k2,-I]] like 05_abelian_extract_k1). *)
sum = FCE[a12 + a1113 + a510];
polk1=SPD[k1,Polarization[k2,-I]];
polk3=SPD[k3,Polarization[k2,-I]];
abK1 = Coefficient[sum /. {polk3->0}, polk1];
Print["MASSIVE abelian k1-coeff LC: ", LeafCount[abK1]];
Print["mq present in abK1? ", !FreeQ[abK1, mq]];
Print["B0 args: ", Union[Cases[abK1,B0[x__]:>{x},Infinity]]];
Print["C0 args: ", Union[Cases[abK1,C0[x__]:>{x},Infinity]]];
Print["D0 args: ", Union[Cases[abK1,D0[x__]:>{x},Infinity]]];
DumpSave["/tmp/ab_coeffK1_massive.mx", abK1];
Print["saved /tmp/ab_coeffK1_massive.mx"];

(* m->0 sanity: master CALL SET should collapse to the massless one *)
abK1m0 = abK1 /. mq->0;
Print["--- m->0 master call set ---"];
Print["B0 args (m->0): ", Union[Cases[abK1m0,B0[x__]:>{x},Infinity]]];
Print["C0 args (m->0): ", Union[Cases[abK1m0,C0[x__]:>{x},Infinity]]];
Print["D0 args (m->0): ", Union[Cases[abK1m0,D0[x__]:>{x},Infinity]]];
