(* TASK 2: the MASSIVE factorizable hard-vertex correction I_{h,f}(Q^2,t).
   Massless reference (2505.10408):
     I_{h,f} = -g^2/(16pi^2) T_i.T_k [ (2Q^2-t) I3^{2m}(Q^2,t) + I2(Q^2) - I2(t) ].
   Here: source(Q^2) -> squark(p1, OFF-shell virtuality t) + antisquark(p2, on-shell m^2),
   one-loop gluon exchange between the two outgoing lines.  Reduce to the validated
   master basis {C0(m^2,t,Q^2; m^2,0,m^2), B0(Q^2;m,m), B0(t;0,m), ...}. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
FCClearScalarProducts[];
ScalarProduct[p1,p1]=tv;            (* emitting leg OFF-shell *)
ScalarProduct[p2,p2]=mq^2;          (* spectator leg on-shell *)
ScalarProduct[p1,p2]=(Qs-tv-mq^2)/2;
dH=InsertFields[CreateTopologies[1,1->2,ExcludeTopologies->{Tadpoles,WFCorrections}],
  {S[1,{ck,cj}]}->{S[13,{1,ci}],-S[13,{1,cl}]},InsertionLevel->{Classes},
  Model->modelM,GenericModel->modelM,ExcludeParticles->{S[1],S[14],F[_]}];
namps=Length[CreateFeynAmp[dH]];
Print["hard-vertex one-loop diagrams: ",namps];
hv=FCFAConvert[CreateFeynAmp[dH],IncomingMomenta->{P},OutgoingMomenta->{p1,p2},
  UndoChiralSplittings->True,ChangeDimension->D,List->False,LoopMomenta->{l},
  SMP->True,Contract->True,DropSumOver->True];
hv=TID[FCE[hv],l,UsePaVeBasis->True,ToPaVe->True];
hv=PaVeReduce[FeynAmpDenominatorExplicit[hv],PaVeAutoReduce->True];
hv=SUNSimplify[hv,SUNNToCACF->True]//Simplify;
Print["masters: ",Union@Cases[hv,(A0|B0|C0|D0)[___],Infinity]];
Print[""];
Print["I_{h,f}^massive (raw, colour kept) ="];
Print[InputForm[Collect[hv/.SMP["g_s"]->gs,{C0[__],B0[__],A0[__]},Simplify]]];
DumpSave["/tmp/ihf_massive.mx",hv];
Print["IHF_DONE"];
