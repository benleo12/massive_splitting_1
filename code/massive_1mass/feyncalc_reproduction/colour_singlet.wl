(* colour_singlet.wl -- shared, CORRECT colour treatment for the radiator amplitudes.

   The source S[1,{ck,cj}] is bi-fundamental, so the amplitude with ck,cj left OPEN is
   NOT gauge invariant: substituting eps^mu -> k2^mu leaves
        T^a_{ci cj} delta_{ck cl}  -  delta_{ci cj} T^a_{ck cl}
   which is precisely the colour-conservation operator T_1^a + T_2^a.  It annihilates
   colour-SINGLET states, so the physical (gauge-invariant) amplitude is the singlet
   projection
        M_singlet = (1/Sqrt[N]) Sum_c M(cj = ck = c) .
   Verified: with this projection k2.M = 0 for BOTH the quark and the squark tree, at
   three independent kinematic points.  With the source open it is nonzero.

   Therefore  |M|^2  and the interference must be built as
        (1/N) [ Sum_c M(c) ] [ Sum_c' M*(c') ]
   with TWO DISTINCT dummies c, c' -- using one dummy for both would wrongly identify
   the source colour of the amplitude and its conjugate.

   NOTE ON A PITFALL: FeynCalc collapses delta_{ab} delta_{ab} -> delta_{ab} rather
   than -> N, so any residual free SUNFDelta after SUNSimplify is a trace waiting to
   be taken and must be replaced by N.  And N[] must NEVER be applied before the
   colour contraction: it turns dummy index labels FCGV[c][1] into FCGV[c][1.] and
   silently destroys index matching.                                              *)

singletM[x_]     := x /. {cj -> csrcA, ck -> csrcA};   (* amplitude        *)
singletMbar[x_]  := x /. {cj -> csrcB, ck -> csrcB};   (* its conjugate    *)

colourContract[x_] := Module[{y},
  y = x /. IndexDelta[a_, b_] :> SUNFDelta[SUNFIndex[a], SUNFIndex[b]];
  y = SUNSimplify[y, Explicit -> True, SUNNToCACF -> False];
  y = SUNSimplify[y, Explicit -> True, SUNNToCACF -> False];
  y = y /. SUNFDelta[_, _] -> SUNN;          (* residual delta = trace -> N *)
  Expand[y /. SUNN -> 3]];

(* build  (1/N) [singlet M] [singlet M*]  and contract all colour *)
singletSquare[M_, Mbar_] := colourContract[ singletM[M] singletMbar[Mbar] ] / 3;

(* gauge check: eps -> k2 must annihilate the singlet-projected amplitude *)
wardSub = Momentum[Polarization[k2, -I], dd_] :> Momentum[k2, dd];
wardResidual[M_] := Simplify[ colourContract[
   DiracSimplify[ ScalarProductExpand[ singletM[M /. wardSub] ] ] ] ];
Print["[colour_singlet] loaded: singletM, singletMbar, colourContract, singletSquare, wardResidual"];
