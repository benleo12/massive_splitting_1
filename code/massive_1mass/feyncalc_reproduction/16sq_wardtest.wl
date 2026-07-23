(* 16sq: DECISIVE gauge test for dropping the closed-squark-loop diagrams.
   15sq found that the 29-diagram squark amplitude VIOLATES the Ward identity
       ((t-mq^2)/2) A + ((u-mq^2)/2) C = 0 .
   Either (a) the 6 dropped closed-squark-loop diagrams {16,17,18,19,30,31} are
   required for gauge invariance -- in which case they cannot be dropped and the
   "apples-to-apples" prescription is gauge-inconsistent -- or (b) the violation
   comes from somewhere else.
   Test: reduce ONLY the 6 dropped diagrams and compute THEIR Ward residual.
     * if  residual(6) = -residual(29)  -> the 6 are required (case a);
     * if  residual(6) = 0              -> the 6 are separately gauge invariant and
                                           the problem is elsewhere (case b).      *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
Get["/tmp/sq_amps_massive.mx"];             (* `amps`, all 35 *)
Print["loaded ", Length[amps], " diagrams"];

upMass[g_] := If[g === 1, mq, 0];
resolveFlavour[e_] := If[FreeQ[e, Index[Generation, 5]], e,
   Sum[ e /. { MQU[Index[Generation,5], _] :> upMass[g],
               MQD[Index[Generation,5], _] :> 0 }, {g, 1, 3}]];
step[a_] := TID[FCE[DiracSimplify[resolveFlavour[a]]], l, UsePaVeBasis->True, ToPaVe->True];
red[a_] := PaVeReduce[FeynAmpDenominatorExplicit[a], PaVeAutoReduce->True];

drop = {16, 17, 18, 19, 30, 31};
totDrop = 0;
Do[Module[{r}, r = Quiet[red[step[amps[[i]]]]]; totDrop = totDrop + r;
   Print["  dropped-diag ", i, " reduced, LC=", LeafCount[r]]], {i, drop}];

kin = {s -> -3, t -> -5, u -> -2, mq -> 1};
M = ScalarProductExpand[FeynAmpDenominatorExplicit[totDrop]] /. kin;
e1 = Pair[Momentum[k1,D],Momentum[Polarization[k2,-I],D]];
e3 = Pair[Momentum[k3,D],Momentum[Polarization[k2,-I],D]];
wardDrop = Simplify[Expand[ScalarProductExpand[
   ((-5-1)/2) Coefficient[M, e1] + ((-2-1)/2) Coefficient[M, e3]] /. kin]];
Print["\nWard residual of the 6 DROPPED squark-loop diagrams:"];
Print[wardDrop];
Print["\nis it zero? ", PossibleZeroQ[wardDrop] || Simplify[wardDrop] === 0];
DumpSave["/tmp/sq_dropped_ward.mx", wardDrop];

(* and the residual of the kept 29, for direct comparison *)
Get["/tmp/sq_reduced_massive.mx"];          (* `tot` = the 29 kept *)
M29 = ScalarProductExpand[FeynAmpDenominatorExplicit[tot]] /. kin;
ward29 = Simplify[Expand[ScalarProductExpand[
   ((-5-1)/2) Coefficient[M29, e1] + ((-2-1)/2) Coefficient[M29, e3]] /. kin]];
Print["\nSUM of the two Ward residuals (must vanish if all 35 are gauge invariant):"];
Print[Simplify[ward29 + wardDrop]];
Print["\nsum is zero? ", PossibleZeroQ[ward29 + wardDrop] ||
                          Simplify[ward29 + wardDrop] === 0];
