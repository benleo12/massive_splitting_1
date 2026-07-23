(* 09ferm (Deliverable B, step 2): resolve the closed-quark-loop FLAVOUR SUM, then
   Dirac-simplify + TID + PaVe-reduce the massive QUARK radiator (19 pure-QCD diagrams).
   INCLUDES the closed quark loops (T_R n_f vacuum-polarisation insertions).

   Flavour sum: FCFAConvert was called with DropSumOver->True, which removes the explicit
   SumOver but leaves the loop-flavour index Index[Generation,5] dangling inside the loop
   masses MQU[Index[Generation,5],..] / MQD[..].  We restore it by summing the diagram over
   the three generations with the model's mass assignments:
        MQU -> {mq, 0, 0}   (gen-1 up-type is the massive quark)
        MQD -> {0, 0, 0}
   Doing this BEFORE the reduction is much cheaper: the massless-loop terms collapse early
   instead of bloating the reduction (the unresolved version reached LC 2.47e6).
   Output: ampRaw_quark_1mass.txt + master basis. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
(* NB: DumpSave stores the SYMBOL `amps`; Get restores it and returns Null. *)
Get["/tmp/ferm_amps_massive.mx"];
Print["loaded ", Length[amps], " quark diagrams"];
If[Length[amps] === 0, Print["ERROR: no diagrams loaded - rerun 08ferm_massive.wl"]; Quit[1]];

(* --- flavour sum over the closed-loop generation --- *)
upMass[g_] := If[g === 1, mq, 0];        (* MQU: gen1 massive, gens 2,3 massless *)
resolveFlavour[e_] := If[FreeQ[e, Index[Generation, 5]], e,
   Sum[ e /. { MQU[Index[Generation,5], _] :> upMass[g],
               MQD[Index[Generation,5], _] :> 0 }, {g, 1, 3}]];

polk2 = Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];
step[a_] := Module[{x},
   x = resolveFlavour[a];
   x = DiracSimplify[x];
   x = TID[FCE[x], l, UsePaVeBasis->True, ToPaVe->True];
   x /. {polk2 -> 0}];
red[a_] := PaVeReduce[FeynAmpDenominatorExplicit[a], PaVeAutoReduce->True];

tot = 0;
Do[Module[{r, hadFl},
   hadFl = !FreeQ[amps[[i]], Index[Generation, 5]];
   r = Quiet[red[step[amps[[i]]]]];
   tot = tot + r;
   Print["  diag ", i, If[hadFl, " (quark loop)", ""], " done, LC=", LeafCount[r]];
 ], {i, Length[amps]}];
Print["TOTAL quark radiator LC: ", LeafCount[tot]];
Print["mq present? ", !FreeQ[tot, mq], "   spinors? ", !FreeQ[tot, Spinor]];
Print["unresolved flavour index left? ", !FreeQ[tot, Index[Generation, 5]]];
Print["A0: ", Union[Cases[tot,A0[x__]:>{x},Infinity]]];
Print["B0: ", Union[Cases[tot,B0[x__]:>{x},Infinity]]];
Print["C0: ", Union[Cases[tot,C0[x__]:>{x},Infinity]]];
Print["D0: ", Union[Cases[tot,D0[x__]:>{x},Infinity]]];
DumpSave["/tmp/ferm_reduced_massive.mx", tot];
Export[FileNameJoin[{base,"for_stefan_1mass","ampRaw_quark_1mass.txt"}],
   ToString[FCE[tot], InputForm], "Text"];
Print["saved ampRaw_quark_1mass.txt"];
