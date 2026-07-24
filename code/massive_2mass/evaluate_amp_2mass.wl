(* ============================================================================
   evaluate_amp_2mass.wl -- driver for the numerical 2-mass squark amplitude.

   Usage:
     AMP_MASTERS=validation/masters.m  wolframscript -file evaluate_amp_2mass.wl
   where masters.m is produced by
     validation/run_masters_2mass.sh <s> <t> <u> <m1sq> <m2sq>  masters.m
   The point (s,t,u,m1sq,m2sq) is embedded in masters.m (key "POINT").
   With no AMP_MASTERS set, uses the bundled validation/masters_demo_2mass.m.
============================================================================ *)
$here = DirectoryName[$InputFileName];
Get[FileNameJoin[{$here, "amplitude_functions_2mass.wl"}]];

mfile = Environment["AMP_MASTERS"];
If[mfile === $Failed || mfile === "",
   mfile = FileNameJoin[{$here, "validation", "masters_demo_2mass.m"}]];

Module[{md, pt, res},
  md = Get[mfile]; pt = md["POINT"];
  If[Head[pt] =!= List,
     Print["ERROR: ", mfile, " has no \"POINT\" -> {s,t,u,m1sq,m2sq}."]; Exit[1]];
  res = amplitudeMass2Ben[md, pt[[1]], pt[[2]], pt[[3]], pt[[4]], pt[[5]]];
  Print["masters: ", mfile];
  Print["amplitudeMass2Ben[s=", pt[[1]], ", t=", pt[[2]], ", u=", pt[[3]],
        ", m1\.b2=", pt[[4]], ", m2\.b2=", pt[[5]], "]  Q2=", pt[[1]]+pt[[2]]+pt[[3]]-pt[[4]]-pt[[5]]];
  Do[Print["  eps^", k-3, " : ", res[[k]]], {k, 1, 6}];
];
