(* Independent AMFlow cross-check of all 13 one-mass masters at a random Euclidean
   point (s=-3, t=-5, u=-2, mq^2=1 -> Q^2=-12).  Auxiliary-mass-flow + DEs + Kira IBP
   (Docker) -- nothing shared with pySecDec.  Convention matches 2F1/pySecDec (ratio=1,
   established on the bubble).  Writes {id -> {eps^kmin..eps^3}} to amflow_1mass_results.m. *)
Get["/Users/user/Library/Wolfram/Applications/AMFlow/AMFlow.m"];
SetReductionOptions["IBPReducer" -> "Kira"];
SetReducerOptions["ReductionMode" -> "FireFly"];
prec = 20; ordEps = 3;

(* point *)
sV=-3; tV=-5; uV=-2; mq=1; QV = sV+tV+uV-2 mq;

(* master specs: {id, legs, conservation, replacement, propagator, target} *)
bub[psq_,ma_,mb_] := {{p1},{},{p1^2->psq},{l^2-ma,(l+p1)^2-mb}, j[fam,1,1]};
tri[P1_,P2_,P3_,ma_,mb_,mc_] := {{k1,k2},{},{k1^2->P1,k2^2->P2,(k1+k2)^2->P3},
   {l^2-ma,(l+k1)^2-mb,(l+k1+k2)^2-mc}, j[fam,1,1,1]};
box[P1_,P2_,P3_,P4_,S12_,S23_,ma_,mb_,mc_,md_] := {{k1,k2,k3},{},
   {k1^2->P1,k2^2->P2,k3^2->P3,(k1+k2)^2->S12,(k2+k3)^2->S23,(k1+k2+k3)^2->P4},
   {l^2-ma,(l+k1)^2-mb,(l+k1+k2)^2-mc,(l+k1+k2+k3)^2-md}, j[fam,1,1,1,1]};

masters = {
 {"B0_t_0_mq",  bub[tV,0,mq]},   {"B0_u_0_mq",  bub[uV,0,mq]},   {"B0_mq2_0_mq", bub[mq,0,mq]},
 {"B0_s_mq_mq", bub[sV,mq,mq]},  {"B0_Q_mq_mq", bub[QV,mq,mq]},
 {"C0_0_mq_t",  tri[0,mq,tV, mq,mq,0]},   {"C0_0_mq_u",  tri[0,mq,uV, mq,mq,0]},
 {"C0_0_s_Q",   tri[0,sV,QV, mq,mq,mq]},  {"C0_mq_mq_s", tri[mq,mq,sV, mq,0,mq]},
 {"C0_mq_t_Q",  tri[mq,tV,QV, mq,0,mq]},  {"C0_mq_u_Q",  tri[mq,uV,QV, mq,0,mq]},
 {"D0_a", box[0,uV,mq,sV, mq,QV, mq,mq,0,mq]},
 {"D0_b", box[mq,tV,QV,sV, 0,mq, mq,0,mq,mq]}};

results = <||>;
Do[Module[{id=m[[1]], spec=m[[2]], legs, cons, repl, prop, tgt, sol, ser, kmin},
   {legs,cons,repl,prop,tgt} = spec;
   ClearAll[AMFlowInfo]; AMFlowInfo["Family"]=fam; AMFlowInfo["Loop"]={l};
   AMFlowInfo["Leg"]=legs; AMFlowInfo["Conservation"]=cons;
   AMFlowInfo["Replacement"]=repl; AMFlowInfo["Propagator"]=prop;
   AMFlowInfo["Numeric"]={}; AMFlowInfo["NThread"]=4;
   Print["=== ", id, " ==="];
   sol = SolveIntegrals[{tgt}, prec, 5];
   ser = tgt /. sol;                              (* Laurent series in eps *)
   results[id] = Table[N[SeriesCoefficient[Series[ser,{eps,0,ordEps}], k], 12], {k, -2, ordEps}];
   Print[id, " -> ", N[results[id], 8]];
 ], {m, masters}];
Put[results, FileNameJoin[{DirectoryName[$InputFileName], "amflow_1mass_results.m"}]];
Print["POINT s=",sV," t=",tV," u=",uV," mq2=",mq," Q2=",QV];
Print["saved amflow_1mass_results.m"];
