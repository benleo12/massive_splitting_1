(* 10ferm (Deliverable B, step 3, FIXED): spin-summed one-loop x tree interference for the
   MASSIVE QUARK radiator, at a NUMERIC kinematic point.

   Fixes over the first attempt (which produced unevaluated garbage):
     * the TREE carries FeynAmpDenominator propagators -> must call
       FeynAmpDenominatorExplicit + ScalarProductExpand so they become explicit
       denominators that the numeric kinematics can evaluate;
     * the Dirac traces generate the dimension D -> substitute D -> 4-2 eps;
     * do NOT take Re[] symbolically (Re/Im of a symbolic expression stays unevaluated) --
       keep the interference linear in the masters and take 2Re only at the very end,
       after the masters have been replaced by numbers (done in 11ferm_evaluate.wl);
     * substitute the numeric kinematics EARLY (collapses 2.1e6 -> ~6e4 leaves).          *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq];
ep = Global`Epsilon;

(* ---------- tree ---------- *)
diagsT = InsertFields[CreateTopologies[0,1->3],
  {S[1,{ck,cj}]}->{F[3,{1,ci}],V[50,{ca}],-F[3,{1,cl}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],S[13],F[10]}];
M0 = Total[FCFAConvert[CreateFeynAmp[diagsT],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->True,LoopMomenta->{},SMP->True,Contract->True,DropSumOver->True]];
M0 = ScalarProductExpand[FeynAmpDenominatorExplicit[M0]];     (* explicit propagators *)
Print["tree LC (explicit props): ", LeafCount[M0]];

(* ---------- one-loop (reduced) ---------- *)
Get["/tmp/ferm_reduced_massive.mx"];      (* restores symbol `tot` *)
M1 = ScalarProductExpand[FeynAmpDenominatorExplicit[tot]];
Print["one-loop LC (symbolic): ", LeafCount[M1]];

(* ---------- numeric kinematics, applied EARLY ---------- *)
kin = {s -> -3, t -> -5, u -> -2, mq -> 1};
M0n = M0 /. kin;  M1n = M1 /. kin;
Print["after numeric kinematics: tree ", LeafCount[M0n], "   1-loop ", LeafCount[M1n]];

(* ---------- WARD IDENTITY CHECK -------------------------------------------
   Write M1^mu eps_mu = A (k1.eps) + B (k2.eps) + C (k3.eps).  Gauge invariance
   requires k2_mu M1^mu = 0, i.e. (k1.k2) A + k2^2 B + (k3.k2) C = 0 with k2^2=0:
        ((t-mq^2)/2) A + ((u-mq^2)/2) C = 0     (as Dirac structures).
   A and C are Dirac-structure valued, so we test the identity after sandwiching
   with the tree and spin-summing.  If it holds, the polarisation sum may be done
   with -g^{mu nu} (no ghosts needed), which is what we do below.            *)
e1 = Pair[Momentum[k1,D],Momentum[Polarization[k2,-I],D]];
e2 = Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];
e3 = Pair[Momentum[k3,D],Momentum[Polarization[k2,-I],D]];
cA = Coefficient[M1n, e1]; cC = Coefficient[M1n, e3];
sand[x_] := Expand[ ScalarProductExpand[ DiracSimplify[
     FermionSpinSum[ x ComplexConjugate[Coefficient[M0n, e1]] ] ] ] /. kin ];
ward = Simplify[ ((-5 - 1)/2) sand[cA] + ((-2 - 1)/2) sand[cC] ];   (* t=-5,u=-2,mq^2=1 *)
Print["WARD  ((t-mq^2)/2) A + ((u-mq^2)/2) C  = ", Simplify[ward]];
Print["WARD identity satisfied? ", TrueQ[Simplify[ward] === 0] || PossibleZeroQ[ward]];

(* ---------- interference: spin + polarisation sums, NO Re[] yet ---------- *)
res = FermionSpinSum[M1n ComplexConjugate[M0n]];
res = DiracSimplify[res];
res = DoPolarizationSums[res, k2, 0];        (* Sum eps^mu eps*^nu -> -g^{mu nu} *)
res = Contract[res];
res = res /. D -> 4 - 2 ep;                  (* dimension from the traces *)
res = ScalarProductExpand[res] /. kin;       (* any kinematics reintroduced by the sums *)
res = Expand[res];
Print["interference LC: ", LeafCount[res]];
Print["leftover FeynAmpDenominator? ", !FreeQ[res, FeynAmpDenominator]];
Print["leftover D?  ", !FreeQ[res, D], "    leftover mq? ", !FreeQ[res, mq]];
Print["leftover Spinor? ", !FreeQ[res, Spinor],
      "    leftover Pair? ", !FreeQ[res, Pair]];
Print["distinct PaVe: ", Length@Union@Cases[res,(A0|B0|C0|D0)[___],Infinity]];
DumpSave["/tmp/quark_interference.mx", res];
Print["saved /tmp/quark_interference.mx"];
