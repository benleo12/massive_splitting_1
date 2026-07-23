(* 14sq (Deliverable B, step 7): drop the closed-SQUARK-loop diagrams, resolve the
   closed-QUARK-loop flavour sum, then TID + PaVe-reduce the SCALAR (squark) radiator.

   Loop content is matched to the quark run (09ferm): gluon + ghost + closed quark
   loops (T_R n_f), no gluinos, and NO closed squark loops.  The discriminator is
   exact and needs no hard-coded indices:
       closed loop of ANY flavour  <=>  the summed index Index[Generation,5] survives
       closed QUARK loop           <=>  additionally a DiracTrace is present
   so  "has Generation-5 index AND no DiracTrace"  is precisely a closed squark loop.
   That drops 6 of the 35 diagrams (16-19, 30, 31), leaving 29.
   Cross-check on the counts: 35 - 27 = 8 quark loops, and 27 - 6 + 8 = 29.        *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
Get["/tmp/sq_amps_massive.mx"];          (* restores the symbol `amps` *)
Print["loaded ", Length[amps], " squark diagrams"];
If[Length[amps] === 0, Print["ERROR: rerun 13sq_generate.wl"]; Quit[1]];

closedSquarkLoop[a_] := !FreeQ[a, Index[Generation, 5]] && FreeQ[a, DiracTrace];
closedQuarkLoop[a_]  := !FreeQ[a, Index[Generation, 5]] && !FreeQ[a, DiracTrace];
drop = Select[Range[Length[amps]], closedSquarkLoop[amps[[#]]] &];
keep = Complement[Range[Length[amps]], drop];
Print["closed SQUARK loops dropped: ", drop, "  (", Length[drop], " diagrams)"];
Print["closed QUARK loops kept    : ",
      Select[Range[Length[amps]], closedQuarkLoop[amps[[#]]] &]];
Print["keeping ", Length[keep], " diagrams"];

(* same flavour assignment as the quark run: MQU gen-1 massive, all else massless *)
upMass[g_] := If[g === 1, mq, 0];
resolveFlavour[e_] := If[FreeQ[e, Index[Generation, 5]], e,
   Sum[ e /. { MQU[Index[Generation,5], _] :> upMass[g],
               MQD[Index[Generation,5], _] :> 0 }, {g, 1, 3}]];

polk2 = Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];
step[a_] := Module[{x},
   x = resolveFlavour[a];
   x = DiracSimplify[x];                       (* traces from the quark loops *)
   x = TID[FCE[x], l, UsePaVeBasis->True, ToPaVe->True];
   x];                                          (* keep eps.k2: needed for the Ward check *)
red[a_] := PaVeReduce[FeynAmpDenominatorExplicit[a], PaVeAutoReduce->True];

tot = 0;
Do[Module[{r},
   r = Quiet[red[step[amps[[i]]]]];
   tot = tot + r;
   Print["  diag ", i, If[closedQuarkLoop[amps[[i]]], " (quark loop)", ""],
         " done, LC=", LeafCount[r]];
 ], {i, keep}];
Print["TOTAL squark radiator LC: ", LeafCount[tot]];
Print["mq present? ", !FreeQ[tot, mq], "   spinors? ", !FreeQ[tot, Spinor]];
Print["unresolved flavour index left? ", !FreeQ[tot, Index[Generation, 5]]];
Print["A0: ", Union[Cases[tot,A0[x__]:>{x},Infinity]]];
Print["B0: ", Union[Cases[tot,B0[x__]:>{x},Infinity]]];
Print["C0: ", Union[Cases[tot,C0[x__]:>{x},Infinity]]];
Print["D0: ", Union[Cases[tot,D0[x__]:>{x},Infinity]]];
DumpSave["/tmp/sq_reduced_massive.mx", tot];
Export[FileNameJoin[{base,"for_stefan_1mass","ampRaw_squark_1mass.txt"}],
   ToString[FCE[tot], InputForm], "Text"];
Print["saved ampRaw_squark_1mass.txt"];
