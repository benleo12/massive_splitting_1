(* 12ferm (Deliverable B, step 5): the spin/polarisation/colour-summed TREE
   |M0|^2 for the massive quark radiator, at the SAME numeric point and with the
   SAME colour machinery as 10/11ferm.  Two purposes:
     (i) diagnose what the leftover SUNFDelta[ci,ck], SUNFDelta[cj,cl] structures
         actually are -- if |M0|^2 carries the same structures, the interference
         normalised by the tree is a clean colour-independent number;
     (ii) provide the denominator for  2Re<M0|M1>/|M0|^2, which is the quantity
          that can be compared between the QUARK and the SCALAR (squark) radiator
          to isolate the magnetic remainder.                                     *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
ep = Global`Epsilon;

diagsT = InsertFields[CreateTopologies[0,1->3],
  {S[1,{ck,cj}]}->{F[3,{1,ci}],V[50,{ca}],-F[3,{1,cl}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],S[13],F[10]}];
M0 = Total[FCFAConvert[CreateFeynAmp[diagsT],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->True,LoopMomenta->{},SMP->True,Contract->True,DropSumOver->True]];
M0 = ScalarProductExpand[FeynAmpDenominatorExplicit[M0]];

Print["--- colour content of the bare tree ---"];
Print[Union@Cases[M0,(SUNT|SUNTF|SUNF|SUNDelta|SUNFDelta|IndexDelta)[___],Infinity]];

kin = {s -> -3, t -> -5, u -> -2, mq -> 1};
M0n = M0 /. kin;

sq = FermionSpinSum[M0n ComplexConjugate[M0n]];
sq = DiracSimplify[sq];
sq = DoPolarizationSums[sq, k2, 0];
sq = Contract[sq];
sq = sq /. D -> 4 - 2 ep;
sq = Expand[ScalarProductExpand[sq] /. kin];

colidx = {ci, cj, ck, cl};
sq = sq /. IndexDelta[a_, b_] /; MemberQ[colidx,a] && MemberQ[colidx,b] :>
             SUNFDelta[SUNFIndex[a], SUNFIndex[b]];
sq = SUNSimplify[sq, Explicit -> True, SUNNToCACF -> False];
sq = SUNSimplify[sq, Explicit -> True, SUNNToCACF -> False];
sq = sq /. SUNFDelta[_, _] -> SUNN;         (* residual delta = trace -> N *)
sq = Expand[sq /. SUNN -> 3];
sq = sq /. SMP["g_s"] -> 1;

Print["\n--- |M0|^2 (spin+pol+colour summed), colour structures that survive ---"];
Print[Union@Cases[sq,(SUNT|SUNTF|SUNF|SUNDelta|SUNFDelta|IndexDelta)[___],Infinity]];
Print["\n|M0|^2 = "];
Print[Simplify[sq]];
Print["\nas a series in eps:"];
Do[Print["  eps^", k, " : ", Simplify[SeriesCoefficient[Series[sq,{ep,0,2}],k]]], {k,0,2}];
DumpSave["/tmp/quark_tree_sq.mx", sq];
Print["\nsaved /tmp/quark_tree_sq.mx"];
