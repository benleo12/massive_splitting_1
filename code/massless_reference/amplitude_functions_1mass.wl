(* ============================================================================
   amplitude_functions_1mass.wl
   ----------------------------------------------------------------------------
   NUMERICAL evaluator for the 1-loop g* -> squark + gluon + antisquark amplitude
   with a MASSIVE squark line (abelian 1/N colour, k1.eps(k2) projection,
   colour/structure-stripped), as an eps-Laurent series {eps^-2 .. eps^3}.

   Design (aligned with the all-numerical / pySecDec strategy):
     * the reduction  ampRaw_1mass.txt  (rational in s,t,u,mq2; VALIDATED via the
       m^2->0 capstone, see validation/capstone_m0.wl) is fixed;
     * the 13 master integrals are supplied NUMERICALLY as bare eps-series
       (convention: bare  Gamma(eps) Int (Symanzik)^{-eps}, identical to pySecDec
       and to masters_bubbles.wl);
     * conversion bare -> the cG/i/mu2 amplitude convention is the loop measure
       factor  i  (proved analytically for the bubble: Bf = i*B0zm; uniform across
       B0/C0/D0), then the global MSBarFac/cGRule is applied at assembly.

   Master numerics:
     bubbles (5)            -> masters_bubbles.wl  (2F1, exact, 16-digit validated)
     triangles+boxes (8)    -> pySecDec  (validation/psd_master.py), per point

   Entry:
     amplitudeMass1Ben[masterBareAssoc, s, t, u, mq2] -> {eps^-2 .. eps^3}
       masterBareAssoc: <| id -> {bare eps-series, from eps^kmin..eps^3} |>
       for all 13 ids (5 bubbles auto-filled from 2F1 if omitted).

   Validation: validation/capstone_m0.wl (m->0 == amplitudeMass0Stefan, eps^3,
   5 points); masters individually vs pySecDec (validation/master_refs.md).
============================================================================ *)

$f1m = DirectoryName[$InputFileName];
ep = Global`Epsilon; ORD = 3;
cGsym = cGsym;  (* symbol *)
cGRule1m = cGsym -> 1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
MSBarFac1m = (16 Pi^2/I) (Exp[EulerGamma]/(4 Pi))^ep;
Imeas = I;   (* bare -> cG/i convention : the loop measure factor i *)

Get[FileNameJoin[{$f1m, "masters_bubbles.wl"}]];   (* bubbleMasters[s,t,u,mq2] -> bare series *)

(* ampRaw_1mass.txt: colour-stripped, in PaVe B0/C0/D0[args] *)
ampRaw1m = ToExpression[Import[FileNameJoin[{$f1m, "ampRaw_1mass.txt"}], "Text"], InputForm] /. {
   SUNN->3, SDF[__]->1, SUNTF[__]->1, SMP["g_s"]->1, ca->1, ci->1, cj->1, ck->1, cl->1};

(* the 13 PaVe -> master-id map (symbolic args; matched BEFORE s,t,u,mq substitution) *)
paveToId = {
  B0[s,mq^2,mq^2]                                  -> "B0_s_mq_mq",
  B0[t,0,mq^2]                                      -> "B0_t_0_mq",
  B0[u,0,mq^2]                                      -> "B0_u_0_mq",
  B0[-2 mq^2+s+t+u,mq^2,mq^2]                       -> "B0_Q_mq_mq",
  B0[mq^2,0,mq^2]                                   -> "B0_mq2_0_mq",
  C0[0,mq^2,t,mq^2,mq^2,0]                          -> "C0_0_mq_t",
  C0[0,mq^2,u,mq^2,mq^2,0]                          -> "C0_0_mq_u",
  C0[mq^2,mq^2,s,mq^2,0,mq^2]                       -> "C0_mq_mq_s",
  C0[0,s,-2 mq^2+s+t+u,mq^2,mq^2,mq^2]              -> "C0_0_s_Q",
  C0[mq^2,t,-2 mq^2+s+t+u,mq^2,0,mq^2]              -> "C0_mq_t_Q",
  C0[mq^2,u,-2 mq^2+s+t+u,mq^2,0,mq^2]              -> "C0_mq_u_Q",
  D0[0,u,mq^2,s,mq^2,-2 mq^2+s+t+u,mq^2,mq^2,0,mq^2]            -> "D0_a",
  D0[mq^2,t,-2 mq^2+s+t+u,s,0,mq^2,mq^2,0,mq^2,mq^2]           -> "D0_b"};

(* bare eps-series, master data entry {kmin, {a_kmin .. a_3}}  (parse_psd.py format) *)
seriesOf[entry_] := Module[{kmin = entry[[1]], lst = entry[[2]]},
  Sum[lst[[k - kmin + 1]] ep^k, {k, kmin, ORD}]];

(* main:  masterData = <| id -> {kmin, {bare eps-series}} |>  for all 13 masters
   (e.g. parsed from pySecDec via validation/parse_psd.py). *)
amplitudeMass1Ben[masterData_, sv_, tv_, uv_, mq2_] := Module[{rules, a},
  rules = (#[[1]] :> Imeas seriesOf[masterData[#[[2]]]]) & /@ paveToId;
  a = ampRaw1m /. rules;                        (* masters -> i * bare series *)
  a = a /. mq -> Sqrt[mq2];                      (* mq^2 -> mq2 in coefficients *)
  a = a /. {s -> sv, t -> tv, u -> uv};
  a = a MSBarFac1m /. cGRule1m;
  Table[N[SeriesCoefficient[Series[a, {ep, 0, ORD}], k], 10], {k, -2, ORD}]];

(* convenience: read parsed pySecDec masters from a .m file and evaluate *)
amplitudeMass1BenFromFile[mfile_, sv_, tv_, uv_, mq2_] :=
  amplitudeMass1Ben[Get[mfile], sv, tv, uv, mq2];

Print["[amplitude_functions_1mass] ready:  amplitudeMass1Ben[masterData, s,t,u,mq2] -> {eps^-2..eps^3}"];
