(* Independent AMFlow cross-check of all 15 TWO-mass masters at a random Euclidean
   point (s=-3, t=-5, u=-2, m1^2=1, m2^2=2 -> Q^2=-13).  Aux-mass-flow + Kira IBP
   (Docker); same bare convention as pySecDec.  Writes amflow_2mass_results.m. *)
Get["/Users/user/Library/Wolfram/Applications/AMFlow/AMFlow.m"];
SetReductionOptions["IBPReducer" -> "Kira"];
SetReducerOptions["ReductionMode" -> "FireFly"];
prec = 20; ordEps = 3;
sV=-3; tV=-5; uV=-2; m1=1; m2=2; QV = sV+tV+uV-m1-m2;

bub[psq_,ma_,mb_] := {{p1},{},{p1^2->psq},{l^2-ma,(l+p1)^2-mb}, j[fam,1,1]};
tri[P1_,P2_,P3_,ma_,mb_,mc_] := {{k1,k2},{},{k1^2->P1,k2^2->P2,(k1+k2)^2->P3},
   {l^2-ma,(l+k1)^2-mb,(l+k1+k2)^2-mc}, j[fam,1,1,1]};
box[P1_,P2_,P3_,P4_,S12_,S23_,ma_,mb_,mc_,md_] := {{k1,k2,k3},{},
   {k1^2->P1,k2^2->P2,k3^2->P3,(k1+k2)^2->S12,(k2+k3)^2->S23,(k1+k2+k3)^2->P4},
   {l^2-ma,(l+k1)^2-mb,(l+k1+k2)^2-mc,(l+k1+k2+k3)^2-md}, j[fam,1,1,1,1]};

masters = {
 {"B0_m1_0_m1", bub[m1,0,m1]},   {"B0_m2_0_m2", bub[m2,0,m2]},
 {"B0_s_m1_m2", bub[sV,m1,m2]},  {"B0_Q_m1_m2", bub[QV,m1,m2]},
 {"B0_t_0_m1",  bub[tV,0,m1]},   {"B0_u_0_m2",  bub[uV,0,m2]},
 {"C0_0_m1_t",  tri[0,m1,tV, m1,m1,0]},   {"C0_0_m2_u",  tri[0,m2,uV, m2,m2,0]},
 {"C0_0_s_Q_a", tri[0,sV,QV, m1,m1,m2]},  {"C0_0_s_Q_b", tri[0,sV,QV, m2,m2,m1]},
 {"C0_m1_m2_s", tri[m1,m2,sV, m1,0,m2]},
 {"C0_m1_u_Q",  tri[m1,uV,QV, m1,0,m2]},  {"C0_m2_t_Q",  tri[m2,tV,QV, m2,0,m1]},
 {"D0_a", box[0,uV,m1,sV, m2,QV, m2,m2,0,m1]},
 {"D0_b", box[m1,tV,QV,sV, 0,m2, m1,0,m1,m2]}};

results = <||>;
Do[Module[{id=m[[1]], spec=m[[2]], legs,cons,repl,prop,tgt, sol, ser},
   {legs,cons,repl,prop,tgt} = spec;
   ClearAll[AMFlowInfo]; AMFlowInfo["Family"]=fam; AMFlowInfo["Loop"]={l};
   AMFlowInfo["Leg"]=legs; AMFlowInfo["Conservation"]=cons;
   AMFlowInfo["Replacement"]=repl; AMFlowInfo["Propagator"]=prop;
   AMFlowInfo["Numeric"]={}; AMFlowInfo["NThread"]=4;
   Print["=== ", id, " ==="];
   sol = SolveIntegrals[{tgt}, prec, 5]; ser = tgt /. sol;
   results[id] = Table[N[SeriesCoefficient[Series[ser,{eps,0,ordEps}], k], 12], {k, -2, ordEps}];
   Print[id, " -> ", N[results[id], 8]];
 ], {m, masters}];
Put[results, FileNameJoin[{DirectoryName[$InputFileName], "amflow_2mass_results.m"}]];
Print["POINT s=",sV," t=",tV," u=",uV," m1sq=",m1," m2sq=",m2," Q2=",QV];
Print["saved amflow_2mass_results.m"];
