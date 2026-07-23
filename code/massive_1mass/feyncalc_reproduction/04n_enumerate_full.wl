(* 04n: enumerate + classify the FULL massive diagram set (abelian + non-abelian)
   for the 1-mass squark process S[1] -> S[13,{1}] V[50] -S[13,{1}].
   Goal (Deliverable A): identify every insertion, its loop content, and whether it
   carries the triple-gluon vertex (SUNF / non-abelian), and locate the external-leg
   self-energy insertions to exclude (amputated, bare MSbar).
   Output: a classification table + the candidate full diagram list. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];

(* full topology set, tadpoles excluded only (we classify WF corrections ourselves) *)
tops = CreateTopologies[1, 1->3, ExcludeTopologies->{Tadpoles}];
diagsl2 = InsertFields[tops,
  {S[1,{cl,cj}]}->{S[13,{1,ci}],V[50,{ca}],S[13,{1,ck}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],F[_]}];
nD = Length[DiagramExtract[diagsl2, All]];
Print["TOTAL insertions: ", nD];

(* also the WFCorrection-free set, to identify self-energy insertions by difference *)
topsNoWF = CreateTopologies[1, 1->3, ExcludeTopologies->{Tadpoles, WFCorrections}];
diagsNoWF = InsertFields[topsNoWF,
  {S[1,{cl,cj}]}->{S[13,{1,ci}],V[50,{ca}],S[13,{1,ck}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],F[_]}];
Print["insertions w/o WF corrections: ", Length[DiagramExtract[diagsNoWF, All]]];

(* classify each insertion of the full set *)
polk2=Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];
conv[idx_]:=FCFAConvert[CreateFeynAmp[DiagramExtract[diagsl2,{idx}]],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->False,LoopMomenta->{l},SMP->True,Contract->False,DropSumOver->True];
Print["idx | LC | SUNF(3g vertex) | #gluonprops | zero?"];
Do[Module[{a, lc, hasF, nglu, zeroQ},
   a = Quiet[conv[i]];
   lc = LeafCount[a];
   hasF = !FreeQ[a, SUNF];
   nglu = Count[a, FeynAmpDenominator[___], Infinity];
   zeroQ = (a === 0);
   Print[i, " | ", lc, " | ", hasF, " | ", nglu, " | ", zeroQ];
 ], {i, nD}];
Print["done"];
