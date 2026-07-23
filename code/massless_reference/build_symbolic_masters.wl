(* build_symbolic_masters.wl — rebuild HF symbolic master Associations from .frm output. *)

SetDirectory[DirectoryName[$InputFileName]];
$baseDir = DirectoryName[$InputFileName];
$hfDir   = FileNameJoin[{ParentDirectory[$baseDir], "MI_via_DE", "hyperform"}];

(* MZV lookup from mzvlow.h. *)
mzvHeaderPath = FileNameJoin[{$hfDir, "src", "mzvlow.h"}];
mzvLines = ReadList[mzvHeaderPath, "String"];
$mzvTable = <||>;
parseMzvRhs[s_String] := ToExpression[StringReplace[s,
  {"z2" -> "Zeta[2]", "z3" -> "Zeta[3]", "z4" -> "Zeta[4]",
   "z5" -> "Zeta[5]", "z6" -> "Zeta[6]", "z7" -> "Zeta[7]",
   "Sinf" -> "0"}]];
Do[
  Module[{m, argStr, rhsStr, argList},
    m = StringCases[ln,
      RegularExpression["Fill\\s+mzv[0-9]+\\s*\\(([0-9,\\s]+)\\)\\s*=\\s*([^;]+)\\s*;"]
        :> {"$1", "$2"}];
    If[Length[m] > 0,
      {argStr, rhsStr} = m[[1]];
      argList = ToExpression["{" <> argStr <> "}"];
      $mzvTable[argList] = parseMzvRhs[rhsStr];
    ];
  ], {ln, mzvLines}];
Print["MZV table: ", Length[$mzvTable], " entries."];

ClearAll[LinfRed];
mzvGet[k_List] := If[KeyExistsQ[$mzvTable, k], $mzvTable[k], Missing["KeyAbsent"]];
LinfRed[args___?(MemberQ[{0, 1}, #] &)] := LinfRed[args] = Module[
  {wl = {args}, n, dual, inputWord, allInserts, count, others, terms, v},
  n = Length[wl];
  Which[
    n == 0, 1,
    AllTrue[wl, # == 0 &] || AllTrue[wl, # == 1 &], 0,
    First[wl] == 0 && Last[wl] == 1,
       v = mzvGet[wl];
       If[!MissingQ[v], v,
         dual = Reverse[1 - wl]; v = mzvGet[dual];
         If[!MissingQ[v], v, Print["MISSING: ", wl]; Linf @@ wl]],
    First[wl] == 0 && Last[wl] == 0, LinfRed @@ Reverse[1 - wl],
    First[wl] == 1,
       inputWord = wl;
       allInserts = Table[Insert[Rest[wl], 1, k], {k, 1, n}];
       count = Count[allInserts, inputWord];
       others = DeleteCases[allInserts, inputWord];
       terms = Total[LinfRed @@@ others];
       If[count == 0, -terms, -terms/count]
  ]];

parseHFOut[file_String, brackSubsts_List] := Module[
   {lines, startIdx, endIdx, block, exprText, expr},
   lines = ReadList[file, "String"];
   startIdx = First @ Flatten @ Position[lines,
     l_String /; StringMatchQ[l,
       RegularExpression["^\\s+[A-Za-z][A-Za-z0-9]*\\s*=\\s*\\(?\\s*$"]], 1, 1];
   endIdx = Last @ Flatten @ Position[StringTrim /@ lines,
     _String?(StringMatchQ[#, RegularExpression["^\\)*\\s*;\\s*$"]] &)];
   block = Take[lines, {startIdx, endIdx}];
   block[[1]] = StringReplace[block[[1]], RegularExpression["^.*=\\s*"] -> ""];
   block[[-1]] = StringReplace[block[[-1]], RegularExpression[";\\s*$"] -> ""];
   exprText = StringJoin[Riffle[block, " "]];
   exprText = StringReplace[exprText,
      {RegularExpression["replace_\\([^)]*\\)"] -> "0",
       RegularExpression["replace_\\[[^]]*\\]"] -> "0"}];
   Do[exprText = StringReplace[exprText, r[[1]] -> r[[2]]], {r, brackSubsts}];
   exprText = StringReplace[exprText, {
      RegularExpression["\\bden\\(([^)]*)\\)"]  :> "den[$1]",
      RegularExpression["\\bLinf\\(([^)]*)\\)"] :> "Linf[$1]",
      RegularExpression["\\bL\\(([^)]*)\\)"]    :> "L[$1]"}];
   expr = ToExpression[exprText, InputForm, Hold];
   ReleaseHold @ (expr /. {
     den[x_] :> 1/x, log[x_] :> Log[x],
     z2 -> Zeta[2], z3 -> Zeta[3], z4 -> Zeta[4],
     z5 -> Zeta[5], z6 -> Zeta[6], z7 -> Zeta[7]})];

c02mPath = FileNameJoin[{$hfDir, "src", "triangle_C0_2m.out"}];
Print["\nReading ", c02mPath, " ..."];
C02m = parseHFOut[c02mPath, {"[t/s]" -> "tau"}];
Print["  parsed; LeafCount = ", LeafCount[C02m]];
C02mCoefs = Association @ Table[k -> Coefficient[C02m, ep, k], {k, -2, 3}];
Print["  per-eps sizes: ", Table[LeafCount[C02mCoefs[k]], {k, -2, 3}]];
Put[C02mCoefs, FileNameJoin[{$baseDir, "C02m_symbolic.m"}]];
Print["  saved C02m_symbolic.m"];

d0Path = FileNameJoin[{$hfDir, "src", "onemass_box_eps3.out"}];
Print["\nReading ", d0Path, " ..."];
D01mass = parseHFOut[d0Path, {"[s/Q2]" -> "sQ", "[u/Q2]" -> "uQ"}];
Print["  parsed; LeafCount = ", LeafCount[D01mass]];
allLinf = Cases[D01mass, _Linf, Infinity] // Union;
D01massRed = D01mass /. (allLinf[[#]] -> (LinfRed @@ allLinf[[#]]) & /@ Range[Length[allLinf]]);
remaining = Cases[D01massRed, _Linf, Infinity] // Union;
Print["  Linf reduced: ", Length[allLinf], " -> ", Length[remaining]];
D01massCoefs = Association @ Table[k -> Coefficient[D01massRed, ep, k], {k, -2, 3}];
Print["  per-eps sizes: ", Table[LeafCount[D01massCoefs[k]], {k, -2, 3}]];
Put[D01massCoefs, FileNameJoin[{$baseDir, "D01mass_symbolic.m"}]];
Print["  saved D01mass_symbolic.m"];
