(* ============================================================================
   amplitude_functions_full.wl  --  assemble the COMPLETE Fig-14 radiator
   (abelian + non-abelian, 18 masters), colour kept symbolic, to extract the
   eps-Laurent pole structure of the full amplitude at a point.
   Entry: fullAmpPoles[masterData, sv,tv,uv,mq2] -> Association eps^k -> colour tensor.
============================================================================ *)
$ff = DirectoryName[$InputFileName];
ep = Global`Epsilon; ORD = 3;
cGsym = cGsym;
cGRuleF = cGsym -> 1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
MSBarFacF = (16 Pi^2/I) (Exp[EulerGamma]/(4 Pi))^ep;
Imeas = I;

ampRawF = ToExpression[Import[FileNameJoin[{$ff, "feyncalc_reproduction",
   "ampRaw_full_1mass.txt"}], "Text"], InputForm] /. {SUNN->3, SMP["g_s"]->1};
(* colour left symbolic: SUNTF, SUNF, SUNFDelta, SDF kept *)

paveToIdF = {
  (* --- 13 abelian --- *)
  B0[s,mq^2,mq^2]->"B0_s_mq_mq", B0[t,0,mq^2]->"B0_t_0_mq", B0[u,0,mq^2]->"B0_u_0_mq",
  B0[-2 mq^2+s+t+u,mq^2,mq^2]->"B0_Q_mq_mq", B0[mq^2,0,mq^2]->"B0_mq2_0_mq",
  C0[0,mq^2,t,mq^2,mq^2,0]->"C0_0_mq_t", C0[0,mq^2,u,mq^2,mq^2,0]->"C0_0_mq_u",
  C0[mq^2,mq^2,s,mq^2,0,mq^2]->"C0_mq_mq_s", C0[0,s,-2 mq^2+s+t+u,mq^2,mq^2,mq^2]->"C0_0_s_Q",
  C0[mq^2,t,-2 mq^2+s+t+u,mq^2,0,mq^2]->"C0_mq_t_Q",
  C0[mq^2,u,-2 mq^2+s+t+u,mq^2,0,mq^2]->"C0_mq_u_Q",
  D0[0,u,mq^2,s,mq^2,-2 mq^2+s+t+u,mq^2,mq^2,0,mq^2]->"D0_a",
  D0[mq^2,t,-2 mq^2+s+t+u,s,0,mq^2,mq^2,0,mq^2,mq^2]->"D0_b",
  (* --- 6 new non-abelian (incl. the A0 tadpole, absent from the abelian sector) --- *)
  A0[mq^2]->"A0_mq",
  C0[0,mq^2,t,0,0,mq^2]->"C0_ggt", C0[0,mq^2,u,0,0,mq^2]->"C0_ggu",
  D0[0,u,-2 mq^2+s+t+u,t,mq^2,mq^2,0,0,mq^2,mq^2]->"D0_c",
  D0[mq^2,0,mq^2,-2 mq^2+s+t+u,t,u,mq^2,0,0,mq^2]->"D0_d",
  D0[mq^2,t,mq^2,u,0,-2 mq^2+s+t+u,0,mq^2,0,mq^2]->"D0_e"};

seriesOf[entry_] := Sum[entry[[2]][[k-entry[[1]]+1]] ep^k, {k, entry[[1]], ORD}];

(* colour-sampling: replace each distinct colour structure by a random value, so the
   pole coefficients become numbers.  If a coefficient vanishes for several independent
   random colour assignments, it vanishes identically. *)
colourSample[expr_, seed_] := Module[{cs, rules},
  cs = Union@Cases[expr, (SUNTF|SUNF|SUNFDelta|SDF)[___], Infinity];
  SeedRandom[seed];
  rules = # -> RandomReal[{-1, 1}, WorkingPrecision -> 30] & /@ cs;
  expr /. rules];

fullAmpPoles[md_, sv_, tv_, uv_, mq2_, nsamp_:3] := Module[{rules, a, ser, coeffs},
  rules = (#[[1]] :> Imeas seriesOf[md[#[[2]]]]) & /@ paveToIdF;
  a = ampRawF /. rules;
  a = a /. mq -> Sqrt[mq2];
  a = a /. {s -> sv, t -> tv, u -> uv};
  a = a MSBarFacF /. cGRuleF;
  ser = Series[a, {ep, 0, ORD}];
  coeffs = Table[SeriesCoefficient[ser, k], {k, -2, ORD}];
  (* for each eps order, evaluate at nsamp random colour assignments *)
  Association @@ Table[
    (k - 3) -> Table[N[colourSample[coeffs[[k]], 100 + j], 12], {j, nsamp}],
    {k, 1, Length[coeffs]}]];
Print["[amplitude_functions_full] ready: fullAmpPoles[masterData, s,t,u,mq2]"];
