(* ===== LiteRed2020 — Evaluation->Quit Kernel FIRST, then run this. =====
   NO Unprotect[j] (it caused j::shdw shadowing). Clean kernel => j lives only in LiteRed`. *)
PrependTo[$Path, "/Users/user/eikonalmasters/litered2020/source"]; << LiteRed2020`

SetDim[d];
Declare[{l, k1, k2}, Vector, {msq, p1s, p2s, ww}, Number];
SetConstraints[{k1, k2}, sp[k1, k1] = p1s; sp[k2, k2] = p2s; sp[k1, k2] = ww];

dir = "/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop/for_stefan_1mass/de_pipeline";
SetDirectory[dir];

NewDsBasis[tri, {sp[l] - msq, sp[l + k1], sp[l + k1 + k2] - msq}, {l},
   Directory -> "tri", FindSymmetries -> True];

(* DIAGNOSTIC: did NewDsBasis build the basis? *)
Print["UniqueSectors: ", UniqueSectors[tri]];
Print["NonZeroSectors: ", NonZeroSectors[tri]];

(* solve, then masters + DE *)
SolvejSector /@ UniqueSectors[tri];
ms = MIs[tri];
Print["Masters: ", ms];  Print["#masters = ", Length[ms]];
de = Outer[IBPReduce[Dinv[#1, #2]] &, ms, {p1s, p2s, ww}, 1];
Export[FileNameJoin[{dir, "c03_DE.m"}], {ms, de}];
Print["Exported c03_DE.m  (", Length[ms], " masters)"];
