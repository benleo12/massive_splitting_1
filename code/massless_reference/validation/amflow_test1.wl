(* AMFlow single-master smoke test: B0(s; mq^2, mq^2) at s=-3, mq^2=1.
   Confirms AMFlow runs here and lets us fix the AMFlow<->pySecDec convention. *)
Get["/Users/user/Library/Wolfram/Applications/AMFlow/AMFlow.m"];
SetReductionOptions["IBPReducer" -> "FIRE+LiteRed"];
AMFlowInfo["Family"] = bmm;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p1^2 -> s};
AMFlowInfo["Propagator"] = {l^2 - msq, (l + p1)^2 - msq};
AMFlowInfo["Numeric"] = {s -> -3, msq -> 1};
AMFlowInfo["NThread"] = 4;
sol = SolveIntegrals[{j[bmm, 1, 1]}, 30, 3];
Print["AMFLOW B0_s_mq_mq = ", sol];
Quit[];
