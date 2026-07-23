(* masters_bubbles.wl -- the 5 massive bubble masters of the 1-mass-squark
   amplitude, as reliable closed forms (2F1 / DE-solution), evaluated as a
   numeric eps-series per kinematic point (same pattern as the massless
   package's hypE box function).  Returns {eps^-1 .. eps^3}.

   VALIDATED against pySecDec to 16 digits at every order, multiple Euclidean
   points (see validation/).  These are classical -> Mathematica-native (no
   handyG needed for the bubbles).

   Normalisation: bare integral  Gamma(eps) Int_0^1 (Symanzik)^{-eps} ...,
   identical convention to pySecDec (no extra factor).  Apply the global
   MSBarFac/c_Gamma at assembly, exactly as in for_stefan/amplitude_functions.wl.
*)

$bubORD = 3;

(* B0(p^2; 0, m^2): one massless + one massive internal line. *)
B0zm[p2_, m2_] := Module[{c = p2/(m2 - p2)},
   Table[N[SeriesCoefficient[Series[
      Gamma[ep] (m2 - p2)^(-ep)/(1 - ep) Hypergeometric2F1[ep, 1 - ep, 2 - ep, -c],
      {ep, 0, $bubORD}], k], 16], {k, -1, $bubORD}]];

(* B0(m^2; 0, m^2): one massless + one massive line, AT p^2 = m^2 (the squark on
   its own mass shell). The 2F1 form B0zm has a threshold pole at p^2=m^2; here the
   Feynman integrand collapses: [m^2(1-x)^2]^{-eps} -> Gamma(eps)(m^2)^{-eps}/(1-2eps). *)
B0zm$onshell[m2_] := Table[N[SeriesCoefficient[Series[
   Gamma[ep] m2^(-ep)/(1 - 2 ep), {ep, 0, $bubORD}], k], 16], {k, -1, $bubORD}];

(* B0(p^2; m^2, m^2): two equal massive internal lines. *)
B0mm[p2_, m2_] := Module[{z = p2/m2},
   Table[N[SeriesCoefficient[Series[
      Gamma[ep] m2^(-ep) (1 - z/4)^(-ep) Hypergeometric2F1[ep, 1/2, 3/2, -z/(4 - z)],
      {ep, 0, $bubORD}], k], 16], {k, -1, $bubORD}]];

(* the five masters of the amplitude, by their (p^2; m1^2,m2^2) signature.
   mq2 = squark mass^2; Q2 = s+t+u-2 mq2. *)
bubbleMasters[s_, t_, u_, mq2_] := With[{Q2 = s + t + u - 2 mq2}, <|
   "B0_t_0_mq"   -> B0zm[t,   mq2],
   "B0_u_0_mq"   -> B0zm[u,   mq2],
   "B0_mq2_0_mq" -> B0zm$onshell[mq2],
   "B0_s_mq_mq"  -> B0mm[s,   mq2],
   "B0_Q_mq_mq"  -> B0mm[Q2,  mq2] |>];
