(* 13sq (Deliverable B, step 6): SCALAR (squark) radiator for the magnetic remainder.
   Process : {S[1,{ck,cj}]} -> {S[13,{1,ci}], V[50,{ca}], -S[13,{1,cl}]}
   Loop content is matched to the QUARK run (08/09ferm, 19 diagrams):
       gluon loop YES, ghost loop YES, closed QUARK loop YES (n_f), gluino NO,
       and closed SQUARK loops MUST BE DROPPED -- they have no counterpart on the
       quark side (there S[13] was excluded), so keeping them would leave a
       squark-loop gluon self-energy contaminating  quark - scalar.
   Same topology exclusions as the quark run: {Tadpoles, WFCorrections}.
   This script generates and CLASSIFIES; it prints the internal field content of
   every diagram so the closed-squark-loop subset can be identified explicitly. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];

diagsS = InsertFields[CreateTopologies[1,1->3,ExcludeTopologies->{Tadpoles,WFCorrections}],
  {S[1,{ck,cj}]}->{S[13,{1,ci}],V[50,{ca}],-S[13,{1,cl}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],F[10]}];

amps = FCFAConvert[CreateFeynAmp[diagsS],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->True,LoopMomenta->{l},SMP->True,Contract->True,DropSumOver->True];
nD = Length[amps];
Print["squark-radiator diagrams (before dropping squark loops): ", nD];

(* A closed fermion loop shows up as a Dirac trace (DiracTrace / no external Spinor).
   A closed SQUARK loop is a pure-scalar bubble on the internal gluon: the diagram
   has NO Dirac structure at all and its loop propagators all carry the squark mass.
   Classify each diagram by (i) presence of a Dirac trace, (ii) the masses appearing
   in the propagators that depend on the loop momentum l.                        *)
loopMasses[a_] := Module[{fads, withL},
  fads = Cases[a, FeynAmpDenominator[x__] :> {x}, Infinity];
  withL = Select[Flatten[fads], !FreeQ[#, l] &];
  Union[Cases[withL, PropagatorDenominator[_, m_] :> m, Infinity]]];
Print["\nidx | LC | DiracTrace? | Spinor? | masses on loop lines"];
Do[Print[i, " | ", LeafCount[amps[[i]]], " | ",
      !FreeQ[amps[[i]], DiracTrace], " | ", !FreeQ[amps[[i]], Spinor], " | ",
      loopMasses[amps[[i]]]], {i, nD}];
Print["\nSMP mq symbol is: ", SMP["m_u"], "   MQU-like symbols seen: ",
      Union[Cases[amps, SMP[_], Infinity]]];
DumpSave["/tmp/sq_amps_massive.mx", amps];
Print["saved /tmp/sq_amps_massive.mx"];
