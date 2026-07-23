(* handyG_wrapper.wl — Mathematica wrapper for handyG's geval Fortran GPL evaluator.

   L[a_1, ..., a_{n-1}, x]  ==  G(a_1, ..., a_{n-1}; x)   (handyG convention).
   geval reads one GPL per input line ("a1,a2,...,a_{n-1},x") and prints one
   result line.  Fortran prints exponents as "E" (e.g. -5.5E-005) and complex
   as "re + im I"; we translate both for Mathematica. *)

(* Locate the geval binary.  IMPORTANT: must be the *patched* geval with a
   >=500-char line buffer (stock handyG geval has character(len=20), which
   truncates the box's long weight-5 GPL input lines and gives wrong results).
   Build it with handyG/build_geval.sh -> handyG/geval_big.  Resolution order:
     1. environment variable  HANDYG_GEVAL  (set to your geval_big path)
     2. handyG/geval_big bundled next to this file
     3. "geval_big" on $PATH
     4. the developer's local build (Ben's machine). *)
$gevalCandidates = {
   Environment["HANDYG_GEVAL"],
   FileNameJoin[{DirectoryName[$InputFileName], "handyG", "geval_big"}],
   "/Users/user/bassi/releases/Bin/local/bin/geval_big"};
$gevalBin = SelectFirst[$gevalCandidates,
   StringQ[#] && FileExistsQ[#] &, "geval_big" (* assume on PATH *)];
If[!FileExistsQ[$gevalBin] && $gevalBin =!= "geval_big",
   Print["[handyG] geval not found; tried ", $gevalCandidates]];
Print["[handyG] geval binary: ", $gevalBin];

(* Lnum (Mathematica NIntegrate GPL evaluator) for the handful of GPLs that
   geval drops at high weight.  Lnum[{a_1..a_{n-1}}, x] = G(a_1..a_{n-1}; x). *)
$LnumPath = FileNameJoin[{DirectoryName[$InputFileName], "Lnum.wl"}];
If[FileExistsQ[$LnumPath], Get[$LnumPath],
   Print["[handyG] Lnum.wl not found at ", $LnumPath, " (no NIntegrate fallback)"]];
lnumOf[L[args__]] := Module[{lst = {args}},
   If[NameQ["Lnum"], Lnum[N[Most[lst], 30], N[Last[lst], 30]], $Failed]];

fmtForGeval[L[args__]] := StringRiffle[
   StringReplace[ToString[CForm[N[#, 18]]], {"*^" -> "e", "e" -> "e"}] & /@ {args}, ","];

(* Fortran numeric token -> Mathematica number (E/D exponent -> *^). *)
toNum[tok_String] := ToExpression[StringReplace[tok, {"E" -> "*^", "e" -> "*^", "D" -> "*^"}]];

parseGevalOutput[line_String] := Module[{trimmed, parts},
   trimmed = StringTrim[line];
   If[trimmed == "", Return["EMPTY"]];
   parts = StringSplit[trimmed];
   Which[
     Length[parts] == 1, toNum[parts[[1]]],
     Length[parts] == 3 && Last[parts] == "I",          (* "re -im I" *)
        toNum[parts[[1]]] + toNum[parts[[2]]] I,
     Length[parts] == 4 && Last[parts] == "I",          (* "re + im I" *)
        toNum[parts[[1]]] + (If[parts[[2]] == "+", 1, -1]) toNum[parts[[3]]] I,
     True, Print["[handyG] unparsed: ", trimmed]; $Failed]];

gevalBatch[lcalls_List] := Module[{tmp, proc, raw, vals},
   (* write to a temp file and run "geval file".  NOTE: geval reads the
      filename into a character(len=20) buffer, so the path must be <= 20
      chars -- use a short /tmp name (calls are sequential). *)
   tmp = "/tmp/hgpl.txt";
   Export[tmp, fmtForGeval /@ lcalls, "Lines"];
   $hgN = If[ValueQ[$hgN], $hgN + 1, 1]; CopyFile[tmp, "/tmp/rb_" <> ToString[$hgN] <> ".txt", OverwriteTarget -> True];
   proc = RunProcess[{$gevalBin, tmp}, All];
   Quiet@DeleteFile[tmp];
   If[proc["ExitCode"] != 0, Print["[handyG] exit ", proc["ExitCode"]]; Return[$Failed]];
   raw = Select[StringSplit[proc["StandardOutput"], "\n"], StringTrim[#] != "" &];
   vals = parseGevalOutput /@ raw;
   If[Length[vals] != Length[lcalls],
     Print["[handyG] ALIGNMENT: ", Length[lcalls], " in, ", Length[vals], " out -> per-call fallback"];
     Return[gevalEach[lcalls]]];
   vals];

(* Safe one-per-call fallback: geval per GPL, and for any geval drops/fails,
   use the Lnum NIntegrate evaluator (handles any weight). *)
gevalEach[lcalls_List] := Module[{nFB = 0, res},
   res = Map[Function[lc, Module[{p = RunProcess[{$gevalBin, "-"}, All, fmtForGeval[lc] <> "\n"], v},
      v = parseGevalOutput[First@StringSplit[p["StandardOutput"] <> "\n", "\n"]];
      If[v === "EMPTY" || v === $Failed,
         nFB++; lnumOf[lc],   (* NIntegrate fallback *)
         v]]], lcalls];
   If[nFB > 0, Print["[handyG] ", nFB, " GPL(s) via Lnum NIntegrate fallback"]];
   res];

evalLwithHandyG[expr_] := Module[{e, lcalls, lvals, bad},
   (* trivial GPLs G(;x)=1 (<=1 entry = no letters); keep only genuine
      numeric-argument GPLs (>=2 entries) for the evaluator. *)
   e = expr /. L[a__] /; Length[{a}] <= 1 :> 1;
   lcalls = DeleteDuplicates@Cases[e, L[a__] /; AllTrue[{a}, NumericQ], Infinity];
   If[Length[lcalls] == 0, Return[e]];
   lvals = gevalBatch[lcalls];          (* file-based geval_big; handles all weights/args *)
   If[lvals === $Failed, Return[$Failed]];
   (* repair only genuine non-numeric ELEMENTS (level {1}; MatchQ excludes the
      List head) via Lnum NIntegrate. *)
   bad = Flatten@Position[lvals, x_ /; MatchQ[x, _?(! NumericQ[#] &)], {1}, Heads -> False];
   If[Length[bad] > 0,
     Print["[handyG] ", Length[bad], " GPL(s) repaired via Lnum NIntegrate"];
     Do[lvals[[i]] = lnumOf[lcalls[[i]]], {i, bad}]];
   e /. Thread[lcalls -> lvals]];

Print["[handyG_wrapper] loaded."];
