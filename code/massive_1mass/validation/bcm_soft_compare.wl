(* ============================================================================
   BCM CROSS-CHECK.  Compare the soft limit of our MASSIVE NON-ABELIAN radiator
   against the one-loop massive soft current of Bierenbaum, Czakon & Mitov
   (arXiv:1107.4384), Case 3 (both emitters massive) -- their function gijCase3.

   Kinematics map (our SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq]  =>
   s=(k1+k3)^2, t=(k1+k2)^2, u=(k2+k3)^2), with p_i = k1, p_j = k3, q = k2:
       mi2 = mj2 = mq^2 ,   pipj = (s-2mq^2)/2 ,
       piq = (t-mq^2)/2  ,  pjq = (u-mq^2)/2 .

   Soft factorisation predicts, for the coefficient of (k1.eps),
       R_soft  =  N(eps) * g_ij / piq ,
   with N(eps) a KINEMATICS-INDEPENDENT normalisation.  So the test is that
       ratio  =  R * piq / g_ij
   is one and the same Laurent series at every soft point and every softness w.
============================================================================ *)
$fn = DirectoryName[$InputFileName];
base = FileNameJoin[{ParentDirectory[ParentDirectory[$fn]]}];
ep = Global`Epsilon; ORD = 3;
Get[FileNameJoin[{ParentDirectory[$fn], "amplitude_functions_nab.wl"}]];

(* ---------------- BCM gijCase3 ---------------------------------------- *)
bcm = FileNameJoin[{base, "refs", "bcm_1107.4384",
                    "OneLoopSoftCurrentBierenbaumCzakonMitov.m"}];
Get[bcm];
Print["BCM file loaded; gijCase3 defined? ", ValueQ[gijCase3] || Head[gijCase3]=!=Symbol];

(* Eq.13 of BCM; only needed from O(e^1) on *)
Fc[x1_?NumericQ, x2_?NumericQ] := Fc[x1,x2] =
  NIntegrate[(x2 Log[1-tt] Log[1-(tt x2)/x1])/(1-tt x2), {tt,0,1},
             WorkingPrecision->20, AccuracyGoal->12, Method->"GlobalAdaptive"];

shorthands = {ai -> (mi2 pjq)/(pipj piq), aj -> (mj2 piq)/(pipj pjq),
              RS -> 4 (mi2 pjq - 2 pipj piq),
              QS -> 16 (mj2 piq^2 - 2 pipj piq pjq + mi2 pjq^2),
              x -> Sqrt[(1-v)/(1+v)], v -> Sqrt[1 - mi2 mj2/pipj^2], IPi -> I Pi};

(* gijCase3 is known through e^1 (ALARMgijC3 multiplies e^2).
   IMPORTANT: truncate the e-series SYMBOLICALLY first.  The O(e) term is the only
   place where Fc[] (a numerical integral) appears, so truncating at e^GORD before
   substituting kinematics keeps Fc out of the way when GORD<=0, and keeps the
   whole substitution cheap. *)
GORD = If[ValueQ[Global`GORDset], Global`GORDset, 0];
gijTrunc = Normal[Series[gijCase3 /. ALARMgijC3 -> 0, {e, 0, GORD}]];
(* BCM Eq.(11):  g_ij(Case 3) = a_S^b ( 2 pipj mu^2 / (2 piq * 2 pjq) )^eps * [ .m file ].
   The .m file gives ONLY the bracket, so the (kinematics-dependent!) prefactor has
   to be restored here -- it is what carries the eps*log(w) of the soft scaling.
   a_S^b = alpha_s^b S_eps/(2 pi) is kinematics-independent and so is absorbed into
   the overall normalisation N(eps) that the test determines. *)
gijAt[sv_, tv_, uv_, mq2_] := Module[{kin, g, pref, PIQ, PJQ, PIPJ},
  PIPJ = (sv - 2 mq2)/2; PIQ = (tv - mq2)/2; PJQ = (uv - mq2)/2;
  kin = {mi2 -> mq2, mj2 -> mq2, pipj -> PIPJ, piq -> PIQ, pjq -> PJQ};
  g = gijTrunc //. shorthands;        (* v feeds x, so iterate *)
  g = g /. kin;
  g = g /. e -> ep;
  pref = Exp[ep Log[(2 PIPJ)/(2 PIQ 2 PJQ)]];    (* mu^2 = 1 *)
  Normal[Series[pref g, {ep, 0, GORD}]]];

(* ---------------- our radiator ---------------------------------------- *)
ourAt[md_, sv_, tv_, uv_, mq2_] := Module[{c},
  c = nabAmp[md, sv, tv, uv, mq2];             (* {eps^-2 .. eps^3} *)
  Sum[c[[k+3]] ep^k, {k, -2, ORD}]];

(* ---------------- the scan ------------------------------------------- *)
mq2 = 1;
pts = {{"A", 20, 2.472135954999579, 6.472135954999579},
       {"B", 40, 8.124555320336759, 4.524555320336759},
       {"C", 100, 7.5505102572168215, 12.449489742783178},
       {"D", 8, 1.8284271247461903, 3.8284271247461903}};
ws  = {0.1, 0.05, 0.02, 0.01};
(* N(eps) = -64 Pi^2 I (1 + eps (a + I b) + ...) is the expected universal
   normalisation.  Extract a,b from the eps^1 coefficient of the ratio:
       ratio_1 = -64 Pi^2 I (a + I b) = 64 Pi^2 b - 64 Pi^2 a I .           *)
NLEAD = -64 Pi^2 I;
abOf[r1_] := {-Im[r1]/(64 Pi^2), Re[r1]/(64 Pi^2)};
tally = {};

Print["\n=================================================================="];
Print["ratio = R_nab * piq / gijCase3   (must be the SAME series everywhere)"];
Print["==================================================================\n"];
Do[Module[{tag=P[[1]], sv=P[[2]], dt=P[[3]], du=P[[4]]},
  Do[Module[{tv, uv, f, md, R, g, piq, rat},
     tv = mq2 + dt w;  uv = mq2 + du w;  piq = (tv - mq2)/2;
     f = FileNameJoin[{"/tmp/softscan", "m_"<>tag<>"_w"<>ToString[w]<>".m"}];
     If[!FileExistsQ[f], Print["  [missing] ", f]; Return[]];
     md = Get[f];
     If[Length[DeleteCases[Keys[md],"POINT"]] < 14,
        Print["  [incomplete] ", f, "  ", Length[Keys[md]]-1, "/14"]; Return[]];
     R = ourAt[md, sv, tv, uv, mq2];
     g = gijAt[sv, tv, uv, mq2];
     rat = Normal[Series[R piq/g, {ep, 0, GORD+2}]];
     Print["point ", tag, "  w=", w, "   s=", sv, " t=", N[tv,8], " u=", N[uv,8]];
     Do[Print["    eps^", k, " : ", Chop[N[SeriesCoefficient[Series[rat,{ep,0,GORD+2}],k],8],10^-7]],
        {k, 0, GORD+2}];
     Print["    R_-2 * piq = ",
           Chop[N[SeriesCoefficient[Series[R,{ep,0,0}],-2] piq, 12]],
           "   (must be w- and s-INDEPENDENT: 32 Pi^2 I = ", N[32 Pi^2, 12], " I)"];
     AppendTo[tally, {tag, sv, w,
        abOf[N[SeriesCoefficient[Series[rat,{ep,0,GORD+2}],1], 12]]}];
     Print[""];
   ], {w, ws}];
 ], {P, pts}];

Print["\n=================================================================="];
Print["N(eps) = -64 Pi^2 I ( 1 + eps (a + I b) + ... )   extracted per point"];
Print["  b should be pi/2 = ", N[Pi/2,10], " universally."];
Print["  a is tested against the hypothesis  a = a0 - (1/2) Log[s]."];
Print["=================================================================="];

Do[Module[{tag=r[[1]], sv=r[[2]], w=r[[3]], a=r[[4,1]], b=r[[4,2]]},
   Print["  ", tag, "  s=", sv, "  w=", w,
         "   a=", N[a,8], "   b=", N[b,8],
         "   a+Log[s]/2=", N[a + Log[sv]/2, 8]]], {r, tally}];
Print["done."];
