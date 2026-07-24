(* Unequal-mass factorizable hard-vertex: source(Q^2) -> squark_i(m1, OFF-shell t)
   + antisquark_k(m2, on-shell m2^2), gluon exchange between the legs.
   Model SQCDBGF_2mass: generation 1 carries m1, generation 2 carries m2. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_2mass","SQCDBGF"}];
FCClearScalarProducts[];
ScalarProduct[p1,p1]=tv;              (* emitting leg (m1) off-shell *)
ScalarProduct[p2,p2]=m2^2;            (* spectator (m2) on-shell *)
ScalarProduct[p1,p2]=(Qs-tv-m2^2)/2;
dH=InsertFields[CreateTopologies[1,1->2,ExcludeTopologies->{Tadpoles,WFCorrections}],
  {S[1,{ck,cj}]}->{S[13,{1,ci}],-S[13,{2,cl}]},InsertionLevel->{Classes},
  Model->modelM,GenericModel->modelM,ExcludeParticles->{S[1],S[14],F[_]}];
Print["diagrams: ",Length[CreateFeynAmp[dH]]];
hv=FCFAConvert[CreateFeynAmp[dH],IncomingMomenta->{P},OutgoingMomenta->{p1,p2},
  UndoChiralSplittings->True,ChangeDimension->D,List->False,LoopMomenta->{l},
  SMP->True,Contract->True,DropSumOver->True];
hv=hv/.{MQU[1]->m1,MQU[2]->m2,MQU[Index[Generation,1]]->m1,MQU[Index[Generation,2]]->m2};
hv=TID[FCE[hv],l,UsePaVeBasis->True,ToPaVe->True];
hv=PaVeReduce[FeynAmpDenominatorExplicit[hv],PaVeAutoReduce->True];
hv=SUNSimplify[hv,SUNNToCACF->True]//Simplify;
cc=SUNTF[{SUNIndex[gg_]},SUNFIndex[ci],SUNFIndex[cj]] SUNTF[{SUNIndex[gg_]},SUNFIndex[ck],SUNFIndex[cl]];
red=Simplify[Coefficient[Expand[hv],
     SUNTF[{SUNIndex[Glu4]},SUNFIndex[ci],SUNFIndex[cj]] SUNTF[{SUNIndex[Glu4]},SUNFIndex[ck],SUNFIndex[cl]]]];
red=16 Pi^2/SMP["g_s"]^2 red;
Print["I_hf^{m1,m2} / (g^2/16pi^2 Ti.Tk) = "];
Print[InputForm[Collect[red,{C0[__],B0[__],A0[__]},Simplify]]];
Print[""];
Print["m2->m1 check vs equal-mass result:"];
eq=Simplify[red/.m2->m1];
Print[InputForm[Collect[eq,{C0[__],B0[__]},Simplify]]];
Print["IHF2_DONE"];
