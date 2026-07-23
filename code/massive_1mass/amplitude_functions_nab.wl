(* ============================================================================
   amplitude_functions_nab.wl -- assemble the MASSIVE NON-ABELIAN radiator.
   Input  : ampRaw_nab_1mass.txt (k1.eps coeff, colour-stripped by nabprefac -> scalar)
   Masters: 14, all already built + AMFlow-validated (subset of the 19-master set).
   Entry  : nabAmp[masterData, s,t,u,mq2] -> {eps^-2 .. eps^3}
   Validated: m->0 reproduces Hoeche's MyBox with N(eps) = -16 Pi^2 I + ... (07nab).
============================================================================ *)
$fn = DirectoryName[$InputFileName];
ep = Global`Epsilon; ORD = 3;
cGsym = cGsym;
cGRuleN = cGsym -> 1/(4 Pi)^(2-ep) Gamma[1+ep] Gamma[1-ep]^2/Gamma[1-2 ep];
MSBarFacN = (16 Pi^2/I) (Exp[EulerGamma]/(4 Pi))^ep;
Imeas = I;

ampRawN = ToExpression[Import[FileNameJoin[{$fn, "ampRaw_nab_1mass.txt"}], "Text"], InputForm];

paveToIdN = {
  B0[mq^2,0,mq^2]                                   -> "B0_mq2_0_mq",
  B0[-2 mq^2+s+t+u,mq^2,mq^2]                       -> "B0_Q_mq_mq",
  C0[0,mq^2,t,0,0,mq^2]                             -> "C0_ggt",
  C0[0,mq^2,t,mq^2,mq^2,0]                          -> "C0_0_mq_t",
  C0[0,mq^2,u,0,0,mq^2]                             -> "C0_ggu",
  C0[0,mq^2,u,mq^2,mq^2,0]                          -> "C0_0_mq_u",
  C0[0,s,-2 mq^2+s+t+u,mq^2,mq^2,mq^2]              -> "C0_0_s_Q",
  C0[mq^2,t,-2 mq^2+s+t+u,mq^2,0,mq^2]              -> "C0_mq_t_Q",
  C0[mq^2,u,-2 mq^2+s+t+u,mq^2,0,mq^2]              -> "C0_mq_u_Q",
  D0[0,u,mq^2,s,mq^2,-2 mq^2+s+t+u,mq^2,mq^2,0,mq^2]            -> "D0_a",
  D0[0,u,-2 mq^2+s+t+u,t,mq^2,mq^2,0,0,mq^2,mq^2]               -> "D0_c",
  D0[mq^2,0,mq^2,-2 mq^2+s+t+u,t,u,mq^2,0,0,mq^2]               -> "D0_d",
  D0[mq^2,t,mq^2,u,0,-2 mq^2+s+t+u,0,mq^2,0,mq^2]               -> "D0_e",
  D0[mq^2,t,-2 mq^2+s+t+u,s,0,mq^2,mq^2,0,mq^2,mq^2]            -> "D0_b"};

seriesOfN[e_] := Sum[e[[2]][[k-e[[1]]+1]] ep^k, {k, e[[1]], ORD}];

nabAmp[md_, sv_, tv_, uv_, mq2_] := Module[{rules, a},
  rules = (#[[1]] :> Imeas seriesOfN[md[#[[2]]]]) & /@ paveToIdN;
  a = ampRawN /. rules;
  a = a /. mq -> Sqrt[mq2];
  a = a /. {s -> sv, t -> tv, u -> uv};
  a = a MSBarFacN /. cGRuleN;
  Table[N[SeriesCoefficient[Series[a, {ep, 0, ORD}], k], 10], {k, -2, ORD}]];
Print["[amplitude_functions_nab] ready: nabAmp[masterData, s,t,u,mq2] -> {eps^-2..eps^3}"];
