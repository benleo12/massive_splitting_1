(* 08ferm (Deliverable B, step 1): MASSIVE QUARK radiator — the fermion analogue of
   the scalar (squark) radiator, needed for the magnetic/spin remainder.
   Process : {S[1,{ck,cj}]} -> {F[3,{1,ci}], V[50,{ca}], -F[3,{1,cl}]}   (massive quarks,
             F[3] carries Mass->MQU = mq in this model, same as the squark).
   Per the vertex decomposition (Eq. 4 of 2505.10408) the quark-gluon vertex splits into
     scalar S^mu  +  magnetic sigma^{mu nu}  +  (EOM term, vanishes on shell),
   so  quark amplitude = scalar radiator (already computed) + MAGNETIC REMAINDER.
   This script generates + reduces the quark amplitude; the remainder is then the
   difference.  Diagnostic only at this stage: report diagram count, Dirac structure,
   and the master basis, to size the computation before committing to the full decomposition (sizing run). *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];

(* PURE QCD: keep external quarks + gluons only.  The SQCD model would also allow
   gluino(F[10])-squark(S[13]) loops on the quark line; those are NOT part of the QCD
   splitting functions of Ref. 2505.10408, so exclude them (and the scalar source
   partners).  Excluding S[13] also guarantees no double counting with the scalar
   radiator, which IS the squark computation. *)
diagsF = InsertFields[CreateTopologies[1,1->3,ExcludeTopologies->{Tadpoles,WFCorrections}],
  {S[1,{ck,cj}]}->{F[3,{1,ci}],V[50,{ca}],-F[3,{1,cl}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],S[13],F[10]}];

amps = FCFAConvert[CreateFeynAmp[diagsF],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->True,LoopMomenta->{l},SMP->True,Contract->True,DropSumOver->True];
nD = Length[amps];
Print["quark-radiator diagrams: ", nD];
Print["idx | LC | has SUNF(nonab)? | has Spinor? | internal content"];
Do[Print[i," | ",LeafCount[amps[[i]]]," | ",!FreeQ[amps[[i]],SUNF]," | ",
         !FreeQ[amps[[i]],Spinor]," | ",
         If[!FreeQ[amps[[i]],SUNF],"nonabelian","abelian"]], {i,nD}];
Print["Dirac heads present: ", Union[Head/@Cases[amps,(DiracGamma|Spinor|DiracSigma)[___],Infinity]]];
DumpSave["/tmp/ferm_amps_massive.mx", amps];
Print["saved /tmp/ferm_amps_massive.mx"];
