(* ============================================================================
   amplitude_functions_2mass.wl
   ----------------------------------------------------------------------------
   NUMERICAL evaluator for the 1-loop g* -> squark + gluon + antisquark amplitude
   with the two external squarks carrying DIFFERENT masses m1, m2 (abelian 1/N,
   k1.eps(k2) projection, colour/structure-stripped), as {eps^-2 .. eps^3}.
   Constraint:  s+t+u = Q^2 + m1^2 + m2^2.

   Same architecture as the 1-mass package: fixed reduction ampRaw_2mass.txt
   (rational in s,t,u,m1^2,m2^2), 15 master integrals supplied numerically as
   bare eps-series (pySecDec / 2F1), converted to the amplitude convention by the
   loop measure factor i, with the global MSBarFac/cGRule applied at assembly.

   Entry:
     amplitudeMass2Ben[masterData, s, t, u, m1sq, m2sq] -> {eps^-2 .. eps^3}
       masterData = <| id -> {kmin, {bare eps-series}} |> for the 15 ids
       (e.g. parsed from pySecDec via validation/parse_psd.py).

   Validation: validation/capstone_m1.wl  (m2->m1 reproduces amplitudeMass1Ben).
============================================================================ *)

$f2m = DirectoryName[$InputFileName];
ep = Global`Epsilon; ORD = 3;
cGsym = cGsym;
cGRule2m = cGsym -> 1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
MSBarFac2m = (16 Pi^2/I) (Exp[EulerGamma]/(4 Pi))^ep;
Imeas = I;

ampRaw2m = ToExpression[Import[FileNameJoin[{$f2m, "ampRaw_2mass.txt"}], "Text"], InputForm] /. {
   SUNN->3, SDF[__]->1, SUNTF[__]->1, SMP["g_s"]->1, ca->1, ci->1, cj->1, ck->1, cl->1};

(* 15 PaVe -> master-id map (symbolic args; matched BEFORE s,t,u,m1,m2 substitution) *)
paveToId2m = {
  B0[m1^2,0,m1^2]                                  -> "B0_m1_0_m1",
  B0[m2^2,0,m2^2]                                  -> "B0_m2_0_m2",
  B0[s,m1^2,m2^2]                                  -> "B0_s_m1_m2",
  B0[-m1^2-m2^2+s+t+u,m1^2,m2^2]                   -> "B0_Q_m1_m2",
  B0[t,0,m1^2]                                      -> "B0_t_0_m1",
  B0[u,0,m2^2]                                      -> "B0_u_0_m2",
  C0[0,m1^2,t,m1^2,m1^2,0]                          -> "C0_0_m1_t",
  C0[0,m2^2,u,m2^2,m2^2,0]                          -> "C0_0_m2_u",
  C0[0,s,-m1^2-m2^2+s+t+u,m1^2,m1^2,m2^2]           -> "C0_0_s_Q_a",
  C0[0,s,-m1^2-m2^2+s+t+u,m2^2,m2^2,m1^2]           -> "C0_0_s_Q_b",
  C0[m1^2,m2^2,s,m1^2,0,m2^2]                       -> "C0_m1_m2_s",
  C0[m1^2,u,-m1^2-m2^2+s+t+u,m1^2,0,m2^2]           -> "C0_m1_u_Q",
  C0[m2^2,t,-m1^2-m2^2+s+t+u,m2^2,0,m1^2]           -> "C0_m2_t_Q",
  D0[0,u,m1^2,s,m2^2,-m1^2-m2^2+s+t+u,m2^2,m2^2,0,m1^2]          -> "D0_a",
  D0[m1^2,t,-m1^2-m2^2+s+t+u,s,0,m2^2,m1^2,0,m1^2,m2^2]          -> "D0_b"};

seriesOf[entry_] := Module[{kmin = entry[[1]], lst = entry[[2]]},
  Sum[lst[[k - kmin + 1]] ep^k, {k, kmin, ORD}]];

amplitudeMass2Ben[masterData_, sv_, tv_, uv_, m1sq_, m2sq_] := Module[{rules, a},
  rules = (#[[1]] :> Imeas seriesOf[masterData[#[[2]]]]) & /@ paveToId2m;
  a = ampRaw2m /. rules;                         (* masters -> i * bare series *)
  a = a /. {m1 -> Sqrt[m1sq], m2 -> Sqrt[m2sq]}; (* m1^2,m2^2 in coefficients *)
  a = a /. {s -> sv, t -> tv, u -> uv};
  a = a MSBarFac2m /. cGRule2m;
  Table[N[SeriesCoefficient[Series[a, {ep, 0, ORD}], k], 10], {k, -2, ORD}]];

amplitudeMass2BenFromFile[mfile_, sv_, tv_, uv_, m1sq_, m2sq_] :=
  amplitudeMass2Ben[Get[mfile], sv, tv, uv, m1sq, m2sq];

Print["[amplitude_functions_2mass] ready: amplitudeMass2Ben[masterData, s,t,u,m1sq,m2sq] -> {eps^-2..eps^3}"];
