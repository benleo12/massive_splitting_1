(* AMFlow cross-check of the 2 NEW BOX masters (B0(0;m,m) is degenerate for AMFlow -- zero-momentum leg --
   and is validated analytically instead against Gamma(eps)(m^2)^-eps) required by the quark (fermion) radiator,
   at s=-3,t=-5,u=-2,mq^2=1.   B0(0;mq2,mq2) and two boxes. *)
Get["/Users/user/Library/Wolfram/Applications/AMFlow/AMFlow.m"];
SetReductionOptions["IBPReducer" -> "Kira"];
SetReducerOptions["ReductionMode" -> "FireFly"];
prec = 20;
sV=-3; tV=-5; uV=-2; mq=1; QV = sV+tV+uV-2 mq;
bub[psq_,ma_,mb_] := {{p1},{},{p1^2->psq},{l^2-ma,(l+p1)^2-mb}, j[fam,1,1]};
box[P1_,P2_,P3_,P4_,S12_,S23_,ma_,mb_,mc_,md_] := {{k1,k2,k3},{},
   {k1^2->P1,k2^2->P2,k3^2->P3,(k1+k2)^2->S12,(k2+k3)^2->S23,(k1+k2+k3)^2->P4},
   {l^2-ma,(l+k1)^2-mb,(l+k1+k2)^2-mc,(l+k1+k2+k3)^2-md}, j[fam,1,1,1,1]};
masters = {
 {"D0_f", box[mq,mq,0,QV, sV,tV, mq,0,mq,mq]},
 {"D0_g", box[mq,mq,0,QV, sV,uV, mq,0,mq,mq]}};
results = <||>;
Do[Module[{id=m[[1]], spec=m[[2]], legs,cons,repl,prop,tgt, sol, ser},
   {legs,cons,repl,prop,tgt} = spec;
   ClearAll[AMFlowInfo]; AMFlowInfo["Family"]=fam; AMFlowInfo["Loop"]={l};
   AMFlowInfo["Leg"]=legs; AMFlowInfo["Conservation"]=cons;
   AMFlowInfo["Replacement"]=repl; AMFlowInfo["Propagator"]=prop;
   AMFlowInfo["Numeric"]={}; AMFlowInfo["NThread"]=4;
   Print["=== ", id, " ==="];
   sol = SolveIntegrals[{tgt}, prec, 5]; ser = tgt /. sol;
   results[id] = Table[N[SeriesCoefficient[Series[ser,{eps,0,3}], k], 12], {k, -2, 3}];
   Print[id, " -> ", N[results[id], 8]];
 ], {m, masters}];
Put[results, FileNameJoin[{DirectoryName[$InputFileName], "amflow_new3_results.m"}]];
Print["saved amflow_new3_results.m"];
