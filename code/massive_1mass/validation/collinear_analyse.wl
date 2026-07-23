(* Assemble the abelian massive radiator along the collinear scan t -> 0
   (s=-3, u=-2 fixed) and print the Laurent coefficients + successive increments,
   which distinguish a finite (mass-regulated) limit from logarithmic growth.
   usage: SCANDIR=<dir> MQ2=<m^2> wolframscript -file collinear_analyse.wl        *)
$d = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[$d], "amplitude_functions_1mass.wl"}]];
dir = Environment["SCANDIR"]; mq2 = ToExpression[Environment["MQ2"]];
If[dir === $Failed, dir = "/tmp/scan_m1"]; If[!NumberQ[mq2], mq2 = 1];
ts = {-1, -0.3, -0.1, -0.03, -0.01, -0.003};
fmt[x_] := ToString[NumberForm[N[x], {9, 4}]];
Print["t-scan, abelian massive radiator, s=-3, u=-2, mq^2=", mq2];
Print["   t         eps^-2     eps^-1     eps^0      (d eps^-1)"];
prev = Null;
Do[Module[{f, md, r, d},
   f = FileNameJoin[{dir, "m_t" <> ToString[t] <> ".m"}];
   If[!FileExistsQ[f], Print["  ", fmt[t], "   (missing)"],
     md = Get[f];
     r = Re[amplitudeMass1Ben[md, -3, t, -2, mq2]];
     d = If[prev === Null, "", fmt[r[[2]] - prev]];
     prev = r[[2]];
     Print["  ", fmt[t], "   ", fmt[r[[1]]], "   ", fmt[r[[2]]], "   ",
           fmt[r[[3]]], "     ", d]]],
 {t, ts}];
Print["(shrinking increments => finite collinear limit; constant increments => log divergence)"];
