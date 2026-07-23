(* Compare AMFlow vs pySecDec for the 5 NEW non-abelian masters at
   s=-3,t=-5,u=-2,mq^2=1. *)
$dir = DirectoryName[$InputFileName];
amf = Get[FileNameJoin[{$dir, "amflow_new5_results.m"}]];   (* id -> {eps^-2..eps^3} *)
psd = Get["/tmp/psd_new5_parsed.m"];                        (* id -> {kmin,{list}} *)
pad[entry_] := Module[{kmin=entry[[1]], lst=entry[[2]]},
  Table[If[k>=kmin, lst[[k-kmin+1]], 0], {k,-2,3}]];
ids = {"C0_ggt","C0_ggu","D0_c","D0_d","D0_e"};
Print["master   max|AMFlow-pySecDec| (eps^-2..eps^3)"];
Print[StringJoin[Table["-",{50}]]];
worst=0;
Do[Module[{a=amf[id], p=pad[psd[id]], d},
   d = Table[Abs[a[[k]]-p[[k]]], {k,1,6}]; worst=Max[worst,Max[d]];
   Print[StringPadRight[id,8],"  ",ScientificForm[N[Max[d]],3],
         "   ",Row[ScientificForm[N[#],2]&/@d," "]]],
 {id, ids}];
Print[StringJoin[Table["-",{50}]]];
Print["WORST over 5 new masters, all orders: ", ScientificForm[N[worst],3]];
