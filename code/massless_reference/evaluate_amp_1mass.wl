(* ============================================================================
   evaluate_amp_1mass.wl -- driver for the numerical 1-mass-squark amplitude.

   Usage:
     wolframscript -file evaluate_amp_1mass.wl <masters.m>

   <masters.m> is the parsed pySecDec master data for a point, produced by
     validation/run_masters.sh  <s> <t> <u> <mq2>  masters.m
   The kinematic point is embedded in <masters.m> (key "POINT" -> {s,t,u,mq2}),
   so it can never mismatch the masters.  Edit the line below to point at a file,
   or pass it via the AMP_MASTERS environment variable.

   Prints amplitudeMass1Ben = {eps^-2 .. eps^3}.
============================================================================ *)
$here = DirectoryName[$InputFileName];
Get[FileNameJoin[{$here, "amplitude_functions_1mass.wl"}]];

mfile = Environment["AMP_MASTERS"];
If[mfile === $Failed || mfile === "",
   mfile = FileNameJoin[{$here, "validation", "masters_demo_asym.m"}]];   (* default: bundled demo *)

Module[{md, pt, res},
  md = Get[mfile];
  pt = md["POINT"];
  If[Head[pt] =!= List,
     Print["ERROR: ", mfile, " has no \"POINT\" -> {s,t,u,mq2} entry."]; Exit[1]];
  res = amplitudeMass1Ben[md, pt[[1]], pt[[2]], pt[[3]], pt[[4]]];
  Print["masters: ", mfile];
  Print["amplitudeMass1Ben[s=", pt[[1]], ", t=", pt[[2]], ", u=", pt[[3]], ", mq2=", pt[[4]],
        "]  Q2=", pt[[1]]+pt[[2]]+pt[[3]]-2 pt[[4]]];
  Do[Print["  eps^", k-3, " : ", res[[k]]], {k, 1, 6}];
];
