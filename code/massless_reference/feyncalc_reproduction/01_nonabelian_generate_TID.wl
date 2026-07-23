(* Step 1: reproduce the 3 non-abelian amp groups through TID (amp[1]). Cache. *)
PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
collab="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop/collaborator_massless";
SetDirectory[collab]; modelAbs=FileNameJoin[{collab,"SQCDBGF","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,0,0];
diags2=InsertFields[CreateTopologies[1,1->3,ExcludeTopologies->{Tadpoles}],
  {S[1,{ck,cj}]}->{S[13,{1,ci}],V[50,{ca}],-S[13,{1,cl}]}, InsertionLevel->{Classes},
  Model->modelAbs, GenericModel->modelAbs, ExcludeParticles->{S[1],S[14],F[_]}];
conv[idx_,lm_]:=FCFAConvert[CreateFeynAmp[DiagramExtract[diags2,{idx}]],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->False,LoopMomenta->lm,SMP->True,Contract->True,
  DropSumOver->True,FinalSubstitutions->{MQU[__]->0}];
tid[a_]:=TID[FCE[a],l,UsePaVeBasis->True,ToPaVe->True]/.{Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]]->0};

(* Cell 26: amp12 (diag 14) *)
amp12$0 = conv[14,{l}];
(* Cell 29: amps1113 (diag 13 + diag 15) with color projections *)
amps1113$0 = Simplify[ScalarProductExpand[
  (conv[13,{-l+k1+k3}] /. {SUNTF[{a_},i_,j_] SUNTF[{a_},m_,l_] SUNTF[{b_},k_,m_] :>
       (I/2) SUNF[b,a,SUNIndex[Glu6]] SUNTF[{SUNIndex[Glu6]},k,l] SUNTF[{a},i,j]})
 +(conv[15,{l+k2}] /. {SUNTF[{a_},i_,j_] SUNTF[{b_},m_,l_] SUNTF[{a_},k_,m_] :>
       (I/2) SUNF[SUNIndex[Glu6],b,a] SUNTF[{a},k,l] SUNTF[{SUNIndex[Glu6]},i,j]})
  ] /. {Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]]->0}];
(* Cell 31: amps510 (diag 5 + diag 12) *)
amps510$0 = {
  SUNSimplify[conv[5,{l}]] /. {SUNTF[{a_,b_},i_,j_] :>
       (I/2) SUNF[a,b,SUNIndex[Glu6]] SUNTF[{SUNIndex[Glu6]},i,j]} /. {Glu5->Glu6,Glu6->Glu5},
  SUNSimplify[conv[12,{l}]] /. {SUNTF[{a_,b_},i_,j_] :>
       (I/2) SUNF[a,b,SUNIndex[Glu6]] SUNTF[{SUNIndex[Glu6]},i,j]} };

Print["amp12$0 LC: ", LeafCount[amp12$0]];
Print["amps1113$0 LC: ", LeafCount[amps1113$0]];
Print["amps510$0 LC: ", LeafCount[Total[amps510$0]]];

amp12$1     = tid[amp12$0];
amps1113$1  = tid[amps1113$0];
amps510$1   = tid[Total[amps510$0]];
Print["TID done. PaVe heads: ",
  Union[Head/@Cases[{amp12$1,amps1113$1,amps510$1},(A0|B0|C0|D0|PaVe)[___],Infinity]]];
DumpSave["/tmp/nab_amp1.mx", {amp12$1, amps1113$1, amps510$1}];
Print["saved /tmp/nab_amp1.mx"];
