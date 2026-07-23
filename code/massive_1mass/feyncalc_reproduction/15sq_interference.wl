(* 15sq (Deliverable B, step 8): SCALAR (squark) one-loop x tree interference AND the
   colour/pol-summed squark tree, at the same numeric point s=-3,t=-5,u=-2,mq^2=1,
   with exactly the same machinery as 10/11/12ferm so that the QUARK and SCALAR
   numbers are directly subtractable:
     * FeynAmpDenominatorExplicit on the tree (explicit propagators);
     * Ward identity check  ((t-mq^2)/2) A + ((u-mq^2)/2) C = 0  -- which is what
       licenses the -g^{mu nu} polarisation sum with no ghosts;
     * DoPolarizationSums[..., k2, 0];
     * IndexDelta -> SUNFDelta, SUNSimplify[Explicit->True], residual delta -> N;
     * NO N[] before the colour contraction (it corrupts dummy index labels).
   Output: /tmp/sq_interference.mx (symbol `res`) and /tmp/sq_tree_sq.mx (symbol `sqtree`). *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
ep = Global`Epsilon;
kin = {s -> -3, t -> -5, u -> -2, mq -> 1};
colidx = {ci, cj, ck, cl};

colourSum[x_] := Module[{y},
  y = x /. IndexDelta[a_, b_] /; MemberQ[colidx,a] && MemberQ[colidx,b] :>
             SUNFDelta[SUNFIndex[a], SUNFIndex[b]];
  y = SUNSimplify[y, Explicit -> True, SUNNToCACF -> False];
  y = SUNSimplify[y, Explicit -> True, SUNNToCACF -> False];
  y = y /. SUNFDelta[_, _] -> SUNN;            (* residual delta is a trace -> N *)
  Expand[y /. SUNN -> 3] /. SMP["g_s"] -> 1];

(* ---------- squark tree ---------- *)
diagsT = InsertFields[CreateTopologies[0,1->3],
  {S[1,{ck,cj}]}->{S[13,{1,ci}],V[50,{ca}],-S[13,{1,cl}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],F[10]}];
M0 = Total[FCFAConvert[CreateFeynAmp[diagsT],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->True,LoopMomenta->{},SMP->True,Contract->True,DropSumOver->True]];
M0 = ScalarProductExpand[FeynAmpDenominatorExplicit[M0]];
Print["squark tree LC: ", LeafCount[M0]];
Print["tree colour content: ",
      Union@Cases[M0,(SUNT|SUNTF|SUNF|SUNDelta|SUNFDelta|IndexDelta)[___],Infinity]];

Get["/tmp/sq_reduced_massive.mx"];        (* restores `tot` *)
M1 = ScalarProductExpand[FeynAmpDenominatorExplicit[tot]];
Print["squark one-loop LC: ", LeafCount[M1]];

M0n = M0 /. kin;  M1n = M1 /. kin;
Print["after numeric kinematics: tree ", LeafCount[M0n], "  1-loop ", LeafCount[M1n]];

(* ---------- Ward identity ---------- *)
e1 = Pair[Momentum[k1,D],Momentum[Polarization[k2,-I],D]];
e3 = Pair[Momentum[k3,D],Momentum[Polarization[k2,-I],D]];
(* For scalars A and C are plain scalar functions (no Dirac structure), so the Ward
   identity can be tested directly, without sandwiching against the tree. *)
cA = Coefficient[M1n, e1]; cC = Coefficient[M1n, e3];
ward = Simplify[ Expand[ScalarProductExpand[
          ((-5 - 1)/2) cA + ((-2 - 1)/2) cC ] /. kin] ];
Print["WARD  ((t-mq^2)/2) A + ((u-mq^2)/2) C = ", Simplify[ward]];
Print["WARD satisfied? ", PossibleZeroQ[ward] || Simplify[ward] === 0];

(* ---------- interference ---------- *)
res = M1n ComplexConjugate[M0n];
res = DoPolarizationSums[res, k2, 0];
res = Contract[res];
res = res /. D -> 4 - 2 ep;
res = Expand[ScalarProductExpand[res] /. kin];
res = colourSum[res];
Print["\ninterference LC: ", LeafCount[res]];
Print["leftover colour: ",
      Union@Cases[res,(SUNT|SUNTF|SUNF|SUNTrace|SUNDelta|SUNFDelta|IndexDelta)[___],Infinity]];
Print["leftover FeynAmpDenominator/D/mq/Pair: ",
      {!FreeQ[res,FeynAmpDenominator], !FreeQ[res,D], !FreeQ[res,mq], !FreeQ[res,Pair]}];
Print["distinct PaVe: ", Length@Union@Cases[res,(A0|B0|C0|D0)[___],Infinity]];
Print["PaVe list (numeric args):"];
Do[Print["   ", InputForm[q]], {q, Union@Cases[res,(A0|B0|C0|D0)[___],Infinity]}];
DumpSave["/tmp/sq_interference.mx", res];

(* ---------- squark tree squared ---------- *)
sqtree = M0n ComplexConjugate[M0n];
sqtree = DoPolarizationSums[sqtree, k2, 0];
sqtree = Contract[sqtree];
sqtree = sqtree /. D -> 4 - 2 ep;
sqtree = Expand[ScalarProductExpand[sqtree] /. kin];
sqtree = colourSum[sqtree];
Print["\nsquark |M0|^2 = ", Simplify[sqtree]];
DumpSave["/tmp/sq_tree_sq.mx", sqtree];
Print["saved /tmp/sq_interference.mx and /tmp/sq_tree_sq.mx"];
