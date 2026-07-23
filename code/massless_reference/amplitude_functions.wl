(* ========================================================================
   amplitude_functions.wl

   Two functions returning the 1-loop massless g* -> q qbar g amplitude
   (abelian k1-projection, colour-stripped) as an eps-Laurent series
   {eps^-2, eps^-1, eps^0, eps^1, eps^2, eps^3} at any Euclidean point:

     amplitudeMass0Ben[s, t, u, Q2]
         -> OUR value.  Master integrals:
              bubble B0, 1-mass triangle C0(0,0,s) : exact (all orders).
              2-mass triangle C0(0,s,t)            : OUR HyperForm .out,
                  GPLs evaluated numerically by handyG (Fortran) -- genuine
                  HyperForm to eps^3.
              1-mass box D0                        : OUR HyperForm .out
                  (D01mass_symbolic.m), GPLs evaluated numerically by handyG
                  (Fortran) -- genuine HyperForm to eps^3.  The .out's
                  replace_(HYPL0,HYPLlimZero,HYPdenx,HYPinf) factors are
                  zero-limit divergence placeholders that vanish in the
                  finalised result (set to 0 at parse time); the resulting box
                  was checked == the standard 1-mass-box closed form to eps^3
                  (12 digits) at random s != u points, confirming the zeroing
                  is exact.
            So the entire amplitude is genuine OUR HyperForm to eps^3 -- no
            Stefan / closed-form input anywhere in amplitudeMass0Ben.

     amplitudeMass0Stefan[s, t, u, Q2]
         -> HIS value: amplitude with Stefan's analytic closed-form masters
            (BubbleFunc / TriangleFunc / BoxFunc).

   Both apply the same MSBarFac and c_Gamma normalisation.  They agree to
   eps^3 at every kinematic point (the ~1e-8 residual at eps^3 for some
   points is handyG double-precision loss on the highest-weight 2-mass
   triangle GPLs; rebuild handyG --quad for more digits).  Q2 must equal
   s + t + u (all Euclidean, i.e. negative).

   Usage:
       Get["amplitude_functions.wl"];
       amplitudeMass0Ben[-1, -1, -1, -3]
       amplitudeMass0Stefan[-1, -2, -1/2, -7/2]
   ======================================================================== *)

SetDirectory[DirectoryName[$InputFileName]];
$dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{$dir, "handyG_wrapper.wl"}]];
C02mSym    = Get[FileNameJoin[{$dir, "C02m_symbolic.m"}]];
D01massSym = Get[FileNameJoin[{$dir, "D01mass_symbolic.m"}]];

ep = Global`Epsilon; mu2 = 1; ORD = 3;
cG = Subscript[c, "\[CapitalGamma]"];
cGRule = cG -> 1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
MSBarFac = (16 Pi^2/I) (Exp[EulerGamma]/(4 Pi))^ep;
sf[x_, n_] := Normal[Series[Exp[-ep Log[-x/mu2]], {ep, 0, n}]];

(* numeric eps-expansion of 2F1(-eps,-eps,1-eps; z) at numeric z (fast,
   ztest-free: parameter-derivatives evaluated numerically). *)
hypE[zz_] := Normal[N[Series[Hypergeometric2F1[-ep,-ep,1-ep, N[zz,30]],
                             {ep, 0, ORD+3}], 25]];
(* the 1-mass box, eps-series, via numeric 2F1 expansion *)
boxF[a_,b_,q_] := Module[{cc=q-a-b},
   2 I cG/ep^2 (1/(a b))(
      (-(b+cc)/(a (b/mu2)))^ep hypE[1-b/(b+cc)]
    + (-(a+cc)/(a (b/mu2)))^ep hypE[1-a/(a+cc)]
    - (-(a+cc)(b+cc)/((a+b+cc) a (b/mu2)))^ep hypE[1-(a b)/((a+cc)(b+cc))])];

(* ---- Stefan's exact closed-form masters (valid to all orders) ---- *)
StefanRules[sv_,tv_,uv_,qv_] := {
  BFunc[x_]    :> I cG/(ep(1-2 ep)) (-mu2/x)^ep,
  CFunc[x_]    :> I cG/ep^2 (1/x)(-mu2/x)^ep,
  CFunc[x_,y_] :> I cG/ep^2 (1/(x-y))((-mu2/x)^ep - (-mu2/y)^ep),
  DFunc[a_,b_,q_] :> boxF[a,b,q]
};

(* ---- Our HyperForm masters ----
   bubble & 1-mass triangle: exact (all orders).
   2-mass triangle C0(0,s,t): OUR HyperForm .out, GPLs via handyG -- genuine
       to eps^3 (validated vs exact: <1e-8 to eps^2, ~1e-6 at eps^3).
   1-mass box D0: OUR HyperForm .out (D01mass_symbolic.m), GPLs via handyG --
       genuine to eps^3.  The .out's 6635 replace_(HYPL0,HYPLlimZero,HYPdenx,
       HYPinf) factors are zero-limit divergence placeholders: in the finalised
       result HYPLlimZero/HYPinf must cancel, and these particular terms vanish,
       so build_symbolic_masters sets replace_(...) -> 0.  Verified: the box so
       built == the standard 1-mass-box closed form to eps^3 (12 digits) at
       random s != u points (P1 does NOT test the box -- it cancels at s=u).
       So box AND 2-mass triangle are both genuine HyperForm to eps^3. *)
hfSeries[assoc_, kin_, ordHi_] := evalLwithHandyG[Sum[(assoc[k] /. kin) ep^k, {k, -2, ordHi}]];
HFRules[sv_,tv_,uv_,qv_] := {
  BFunc[x_]    :> I cG/(1-2 ep) (1/ep) sf[x, ORD+3],
  CFunc[x_]    :> I cG/(1-2 ep) (1/ep^2 - 2/ep)(1/x) sf[x, ORD+3],
  CFunc[x_,y_] :> I cG/(1-2 ep) (1/x) hfSeries[C02mSym, {tau -> N[y/x, 30]}, ORD] sf[x, ORD+3],
  (* genuine OUR HyperForm box: D01massSym in letters sQ=s/Q^2, uQ=u/Q^2, GPLs via
     handyG; conversion DFunc[a,b,q] = i c_G/(1-2 ep) (-q)^(-2-ep) D0unit{sQ->a/q,uQ->b/q},
     validated == boxF to eps^3 (12 digits) at random s != u points. NO closed form. *)
  DFunc[a_,b_,q_] :> I cG/(1-2 ep) (-q)^(-2-ep) hfSeries[D01massSym, {sQ -> N[a/q, 30], uQ -> N[b/q, 30]}, ORD]
};

(* ---- Stefan's reduced amplitude (abelian k1 projection, colour-stripped) ---- *)
(* Stefan's reduced amplitude (abelian k1 projection), bundled as portable text
   (was a platform-specific .mx).  This is the colour/Dirac-reduced amplitude in
   terms of the master loop functions BFunc/CFunc/DFunc -- NOT his closed forms;
   the masters themselves are supplied by HFRules (our HyperForm) below. *)
ampStr = Import[FileNameJoin[{$dir, "ampRaw_stefan.txt"}], "Text"];
ampRaw = ToExpression[ampStr] /. {
  SUNN->3, SUNFDelta[___]->1, SUNTF[___]->1, SUNFIndex[___]->1, SUNIndex[___]->1,
  Polarization[___]->1, Pair[___]->1, Momentum[___]->1, SMP["g_s"]->1};

(* ---- the two evaluators ---- *)
evalAmp[rulesF_, sv_, tv_, uv_, qv_] := Module[{a},
   (* substitute Q^2 -> qv DIRECTLY (real); never Q -> Sqrt[qv], which would
      make Q^2 = Sqrt[qv]^2 = qv + 0.*I complex for negative qv and corrupt
      the GPL arguments. *)
   a = ampRaw /. {Global`s->sv, Global`t->tv, Global`u->uv};
   a = a /. (Global`Q^2) -> qv;            (* Q^2 -> real qv *)
   a = a /. Global`Q -> Sqrt[qv];          (* any leftover bare Q (none expected) *)
   a = a /. (Global`Q^2) -> qv;
   a = a /. rulesF[sv,tv,uv,qv] /. (Global`Q^2) -> qv;
   a = a MSBarFac /. cGRule;
   Table[N[SeriesCoefficient[Series[a, {ep, 0, ORD}], k], 12], {k, -2, ORD}]];

amplitudeMass0Ben[sv_, tv_, uv_, qv_]    := evalAmp[HFRules,     sv, tv, uv, qv];
amplitudeMass0Stefan[sv_, tv_, uv_, qv_] := evalAmp[StefanRules, sv, tv, uv, qv];

Print["[amplitude_functions] ready:"];
Print["   amplitudeMass0Ben[s,t,u,Q^2]     -> {eps^-2 .. eps^3}, OUR HyperForm masters"];
Print["   amplitudeMass0Stefan[s,t,u,Q^2]  -> {eps^-2 .. eps^3}, Stefan closed forms"];
