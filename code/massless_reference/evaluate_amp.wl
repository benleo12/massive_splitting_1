(* ========================================================================
   evaluate_amp.wl  -- batch/report driver.

   Thin wrapper over amplitude_functions.wl (the canonical definitions of
   amplitudeMass0Ben / amplitudeMass0Stefan).  Prints, at several Euclidean
   points, OUR genuine HyperForm amplitude vs Stefan's analytic closed forms,
   to eps^3, with the per-order difference.

   Usage:
       wolframscript -file evaluate_amp.wl
   Interactive:
       Get["evaluate_amp.wl"];  reportAt[s,t,u,Q2]    (Q2 = s+t+u)
   ======================================================================== *)

SetDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "amplitude_functions.wl"}]];

fmt[c_] := Module[{r = N[Re[c], 9], i = N[Im[c], 9]},
  Which[Abs[r] < 10^-7 && Abs[i] > 10^-7, ToString[CForm[i]] <> " i",
        Abs[i] < 10^-7 && Abs[r] > 10^-7, ToString[CForm[r]],
        Abs[r] < 10^-7 && Abs[i] < 10^-7, "~0",
        True, ToString[CForm[r]] <> "+" <> ToString[CForm[i]] <> "i"]];

reportAt[s_, t_, u_, q2_] := Module[{ben, stef},
  Print[""]; Print[StringRepeat["=", 78]];
  Print["  s=", s, " t=", t, " u=", u, " Q^2=", q2];
  Print[StringRepeat["=", 78]];
  ben  = amplitudeMass0Ben[s, t, u, q2];
  stef = amplitudeMass0Stefan[s, t, u, q2];
  Print["  order   OUR (HyperForm masters)        Stefan (closed forms)         |diff|"];
  Do[Print["  eps^", If[k >= 0, " ", ""], k, "   ",
       StringPadRight[fmt[ben[[k + 3]]], 30], StringPadRight[fmt[stef[[k + 3]]], 30],
       ToString[CForm[N[Abs[ben[[k + 3]] - stef[[k + 3]]], 6]]]], {k, -2, 3}];
];

Print[""]; Print[StringRepeat["#", 78]];
Print["#  Amplitude (abelian k1 projection): OUR HyperForm masters vs Stefan closed forms"];
Print["#  -> eps^-2 .. eps^3, both through the identical MSBarFac and c_Gamma."];
Print[StringRepeat["#", 78]];

reportAt[-1, -1, -1, -3];
reportAt[-2, -1, -1, -4];
reportAt[-1, -2, -1/2, -7/2];
reportAt[-3, -2, -2, -7];

Print[""]; Print[StringRepeat["#", 78]];
Print["#  For your own point:  reportAt[s,t,u,Q^2]   (Q^2 = s+t+u, all negative)"];
Print[StringRepeat["#", 78]];
