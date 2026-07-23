(* AMFlow cross-check of the 5 NEW non-abelian masters at s=-3,t=-5,u=-2,mq^2=1. *)
Get["/Users/user/Library/Wolfram/Applications/AMFlow/AMFlow.m"];
SetReductionOptions["IBPReducer" -> "Kira"];
SetReducerOptions["ReductionMode" -> "FireFly"];
prec = 20;
sV=-3; tV=-5; uV=-2; mq=1; QV = sV+tV+uV-2 mq;
tri[P1_,P2_,P3_,ma_,mb_,mc_] := {{k1,k2},{},{k1^2->P1,k2^2->P2,(k1+k2)^2->P3},
   {l^2-ma,(l+k1)^2-mb,(l+k1+k2)^2-mc}, j[fam,1,1,1]};
box[P1_,P2_,P3_,P4_,S12_,S23_,ma_,mb_,mc_,md_] := {{k1,k2,k3},{},
   {k1^2->P1,k2^2->P2,k3^2->P3,(k1+k2)^2->S12,(k2+k3)^2->S23,(k1+k2+k3)^2->P4},
   {l^2-ma,(l+k1)^2-mb,(l+k1+k2)^2-mc,(l+k1+k2+k3)^2-md}, j[fam,1,1,1,1]};
masters = {
 {"C0_ggt", tri[0,mq,tV, 0,0,mq]}, {"C0_ggu", tri[0,mq,uV, 0,0,mq]},
 {"D0_c", box[0,uV,QV,tV, mq,mq, 0,0,mq,mq]},
 {"D0_d", box[mq,0,mq,QV, tV,uV, mq,0,0,mq]},
 {"D0_e", box[mq,tV,mq,uV, 0,QV, 0,mq,0,mq]}};
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
Put[results, FileNameJoin[{DirectoryName[$InputFileName], "amflow_new5_results.m"}]];
Print["saved amflow_new5_results.m"];
