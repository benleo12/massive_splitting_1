(* Compare AMFlow (independent: aux-mass-flow + Kira IBP) vs pySecDec, all 13
   one-mass masters at s=-3,t=-5,u=-2,mq^2=1.  Both in the same bare convention
   (ratio=1 established on the bubble).  Prints max|AMFlow - pySecDec| per master. *)
$dir = DirectoryName[$InputFileName];
amf = Get[FileNameJoin[{$dir, "amflow_1mass_results.m"}]];   (* id -> {eps^-2..eps^3} *)
psd = Get["/tmp/psd_cmp_point.m"];                            (* id -> {kmin,{list}}, + POINT *)

pad[entry_] := Module[{kmin=entry[[1]], lst=entry[[2]]},     (* pySecDec -> {eps^-2..eps^3} *)
  Table[If[k>=kmin, lst[[k-kmin+1]], 0], {k,-2,3}]];

ids = Sort[DeleteCases[Keys[psd], "POINT"]];
Print["master           max|AMFlow-pySecDec|   (eps^-2 .. eps^3)"];
Print[StringJoin[Table["-",{58}]]];
worst = 0;
Do[Module[{a=amf[id], p=pad[psd[id]], d},
   d = Table[Abs[a[[k]]-p[[k]]], {k,1,6}];
   worst = Max[worst, Max[d]];
   Print[StringPadRight[id,16], "  ", ScientificForm[N[Max[d]],3],
         "   ", Row[ScientificForm[N[#],2]&/@d, " "]];
 ], {id, ids}];
Print[StringJoin[Table["-",{58}]]];
Print["WORST over all 13 masters, all orders: ", ScientificForm[N[worst],3]];
