(* 17redo: recompute BOTH the quark and the squark tree-squared and one-loop x tree
   interference with the CORRECT colour-singlet source projection (see colour_singlet.wl).
   Supersedes the numbers from 11ferm/12ferm/15sq, which used an open bi-fundamental
   source -- a projection under which the amplitude is not even gauge invariant.

   Dropping the eps.k2 terms of the one-loop amplitude (09ferm does this) is harmless
   HERE: it changes the interference by B times (k2 dot M0bar), and that contraction
   vanishes by the tree Ward identity, which the singlet projection now satisfies. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
fcr = FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction"}];
Get[FileNameJoin[{fcr, "colour_singlet.wl"}]];
modelM=FileNameJoin[{fcr,"SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
ep = Global`Epsilon;
kin = {s -> -3, t -> -5, u -> -2, mq -> 1};

tree[out_, excl_] := Module[{d, a},
  d = InsertFields[CreateTopologies[0,1->3], {S[1,{ck,cj}]}->out, InsertionLevel->{Classes},
        Model->modelM, GenericModel->modelM, ExcludeParticles->excl];
  a = Total[FCFAConvert[CreateFeynAmp[d], IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},
        UndoChiralSplittings->True, ChangeDimension->D, List->True, LoopMomenta->{},
        SMP->True, Contract->True, DropSumOver->True]];
  ScalarProductExpand[FeynAmpDenominatorExplicit[a]]];

finish[x_, fermionic_] := Module[{y},
  y = If[fermionic, DiracSimplify[FermionSpinSum[x]], x];
  y = DoPolarizationSums[y, k2, 0];
  y = Contract[y];
  y = y /. D -> 4 - 2 ep;
  y = Expand[ScalarProductExpand[y] /. kin];
  colourContract[y]/3];

doit[lbl_, out_, excl_, loopfile_, fermionic_] := Module[{M0, M0n, M1, M1n, sq, inter},
  Print["\n=========== ", lbl, " ==========="];
  M0 = tree[out, excl];  M0n = M0 /. kin;
  Print["tree Ward residual (singlet, eps->k2): ", wardResidual[M0 /. kin]];
  sq = finish[ singletM[M0n] singletMbar[ComplexConjugate[M0n]], fermionic ];
  Print["|M0|^2 (singlet) = ", Simplify[sq /. SMP["g_s"] -> 1]];
  Get[loopfile];                                  (* restores `tot` *)
  M1 = ScalarProductExpand[FeynAmpDenominatorExplicit[tot]];
  M1n = M1 /. kin;
  inter = finish[ singletM[M1n] singletMbar[ComplexConjugate[M0n]], fermionic ];
  inter = inter /. SMP["g_s"] -> 1;
  Print["interference LC: ", LeafCount[inter]];
  Print["leftover colour: ",
    Union@Cases[inter,(SUNT|SUNTF|SUNF|SUNTrace|SUNDelta|SUNFDelta|IndexDelta)[___],Infinity]];
  Print["distinct PaVe: ", Length@Union@Cases[inter,(A0|B0|C0|D0)[___],Infinity]];
  {sq /. SMP["g_s"] -> 1, inter}];

qk = doit["QUARK", {F[3,{1,ci}],V[50,{ca}],-F[3,{1,cl}]}, {S[1],S[14],S[13],F[10]},
          "/tmp/ferm_reduced_massive.mx", True];
treeQ = qk[[1]]; interQ = qk[[2]];
DumpSave["/tmp/quark_singlet.mx", {treeQ, interQ}];

sqk = doit["SQUARK", {S[13,{1,ci}],V[50,{ca}],-S[13,{1,cl}]}, {S[1],S[14],F[10]},
           "/tmp/sq_reduced_massive.mx", False];
treeS = sqk[[1]]; interS = sqk[[2]];
DumpSave["/tmp/squark_singlet.mx", {treeS, interS}];
Print["\nsaved /tmp/quark_singlet.mx and /tmp/squark_singlet.mx"];
