(* 04n v2: generate the FULL diagram set (physical, amputated: Tadpoles+WFCorrections
   excluded), keep FULL colour (abelian + non-abelian), count + classify, then do the
   k1 AND k3 projections and TID/PaVe reduce -> ampRaw_full_1mass.txt (complete Fig-14
   radiator, both colour structures, both projections).  Masters are the SAME as the
   abelian case (fixed by the loop topology). *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];

tops = CreateTopologies[1, 1->3, ExcludeTopologies->{Tadpoles, WFCorrections}];
diags = InsertFields[tops,
  {S[1,{cl,cj}]}->{S[13,{1,ci}],V[50,{ca}],S[13,{1,ck}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],F[_]}];

amps = FCFAConvert[CreateFeynAmp[diags],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->True,LoopMomenta->{l},SMP->True,Contract->True,DropSumOver->True];
nD = Length[amps];
Print["TOTAL physical diagrams: ", nD];
Print["idx | LC | non-abelian(SUNF)? | nonzero?"];
Do[Print[i," | ",LeafCount[amps[[i]]]," | ",!FreeQ[amps[[i]],SUNF]," | ",amps[[i]]=!=0], {i,nD}];
nNab = Count[amps, a_/;!FreeQ[a,SUNF]];
Print["# non-abelian (triple-gluon-vertex) diagrams: ", nNab];
Print["# abelian-only diagrams: ", nD-nNab];
DumpSave["/tmp/full_amps_1mass.mx", amps];
Print["saved /tmp/full_amps_1mass.mx"];
