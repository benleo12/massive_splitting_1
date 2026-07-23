(* ============================================================================
   CAPSTONE 2 for the 2-mass amplitude (NUMERICAL m2->m1 nesting):
   amplitudeMass2Ben at m1^2 = m2^2 = 1  ==  amplitudeMass1Ben at mq^2 = 1,
   at the SAME asymmetric kinematic point (s=-2,t=-3,u=-1, Q2=-8), where the box
   contributes non-trivially.
   Needs:
     for_stefan_2mass/validation/masters_equalmass.m   (15 masters at m1^2=m2^2=1)
     for_stefan_1mass/validation/masters_demo_asym.m   (13 masters at mq^2=1)
   Run:  wolframscript -file capstone_m1_2mass.wl
============================================================================ *)
$dir = DirectoryName[$InputFileName]; $f2 = ParentDirectory[$dir];
$f1 = FileNameJoin[{ParentDirectory[$f2], "for_stefan_1mass"}];

Get[FileNameJoin[{$f2, "amplitude_functions_2mass.wl"}]];
md2 = Get[FileNameJoin[{$dir, "masters_equalmass.m"}]];
a2  = amplitudeMass2Ben[md2, -2, -3, -1, 1, 1];

Get[FileNameJoin[{$f1, "amplitude_functions_1mass.wl"}]];
md1 = Get[FileNameJoin[{$f1, "validation", "masters_demo_asym.m"}]];
a1  = amplitudeMass1Ben[md1, -2, -3, -1, 1];

Print["m2->m1 nesting at s=-2,t=-3,u=-1, mass^2=1 (Q2=-8):"];
Print["  amplitudeMass2Ben(m1=m2=1) : ", N[a2, 6]];
Print["  amplitudeMass1Ben(mq2=1)   : ", N[a1, 6]];
Print["  |diff| per order eps^-2..3 : ", Table[ScientificForm[N[Abs[a2[[k]]-a1[[k]]]],3], {k,1,6}]];
