(* 11ferm (Deliverable B, step 4): substitute the validated master values into the
   spin-summed quark interference and extract its eps-Laurent series at the point
   s=-3, t=-5, u=-2, mq^2=1  (so Q^2 = -12).
   The interference was built with NUMERIC kinematics, so the PaVe arguments are numbers;
   the map below is the numeric-argument version of the master dictionary. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]]; PrependTo[$Path, "/tmp/fc93"];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0;
ep = Global`Epsilon; ORD = 3;
Get["/tmp/quark_interference.mx"];      (* restores symbol `res` *)
Print["interference LC: ", LeafCount[res]];

md = Join[Get["/tmp/masters19.m"], Get["/tmp/psd_new3_parsed.m"]];
Print["master data entries: ", Length[DeleteCases[Keys[md],"POINT"]]];
ser[id_] := Module[{e = md[id]}, Sum[e[[2]][[k-e[[1]]+1]] ep^k, {k, e[[1]], ORD}]];
I0 = I;   (* bare -> amplitude convention: loop measure factor i *)

num = {
 A0[1] -> "A0_mq",
 B0[0,1,1] -> "B0_0_mq_mq", B0[1,0,1] -> "B0_mq2_0_mq", B0[-3,1,1] -> "B0_s_mq_mq",
 B0[-5,0,1] -> "B0_t_0_mq", B0[-2,0,1] -> "B0_u_0_mq", B0[-12,1,1] -> "B0_Q_mq_mq",
 C0[0,1,-5,0,0,1] -> "C0_ggt",   C0[0,1,-5,1,1,0] -> "C0_0_mq_t",
 C0[0,1,-2,0,0,1] -> "C0_ggu",   C0[0,1,-2,1,1,0] -> "C0_0_mq_u",
 C0[0,-3,-12,1,1,1] -> "C0_0_s_Q", C0[1,1,-3,1,0,1] -> "C0_mq_mq_s",
 C0[1,-5,-12,1,0,1] -> "C0_mq_t_Q", C0[1,-2,-12,1,0,1] -> "C0_mq_u_Q",
 D0[0,-2,1,-3,1,-12,1,1,0,1] -> "D0_a",   D0[0,-2,-12,-5,1,1,0,0,1,1] -> "D0_c",
 D0[1,0,1,-12,-5,-2,1,0,0,1] -> "D0_d",   D0[1,1,0,-12,-3,-5,1,0,1,1] -> "D0_f",
 D0[1,1,0,-12,-3,-2,1,0,1,1] -> "D0_g",   D0[1,-5,1,-2,0,-12,0,1,0,1] -> "D0_e",
 D0[1,-5,-12,-3,0,1,1,0,1,1] -> "D0_b"};

(* report which PaVe actually occur, and whether every one is covered by the map *)
occ = Union[Cases[res, (A0|B0|C0|D0)[___], Infinity]];
Print["distinct PaVe in interference: ", Length[occ]];
unmapped = Select[occ, FreeQ[num[[All,1]], #] &];
Print["UNMAPPED (must be empty): ", unmapped];

(* ---- COLOUR SUM ----------------------------------------------------------
   Sum over all colours: quark ci, antiquark cl, gluon ca, and the two source
   indices ck,cj (each appears once in M1 and once in M0*, so they contract).
   Three things have to be got right:
     * FeynArts emits IndexDelta[..] for the colour deltas of the source vertex;
       FeynCalc's SUN machinery only contracts SUNFDelta, so convert them;
     * SUNTrace is left symbolic unless Explicit->True;
     * NEVER apply N[] before the contraction -- it turns dummy index labels
       FCGV[c][1] into FCGV[c][1.] and silently destroys index matching.      *)
colidx = {ci, cj, ck, cl};
col = res /. IndexDelta[a_, b_] /; MemberQ[colidx, a] && MemberQ[colidx, b] :>
                SUNFDelta[SUNFIndex[a], SUNFIndex[b]];
col = SUNSimplify[col, Explicit -> True, SUNNToCACF -> False];
col = SUNSimplify[col, Explicit -> True, SUNNToCACF -> False];
(* FeynCalc collapses  delta_{ab} delta_{ab} -> delta_{ab}  instead of -> N.
   The residual free delta is therefore a TRACE waiting to be taken: replace it by N.
   Verified against the hand colour algebra of the tree, where
     <C1 C1*> = Tr(T^a T^a) * sum_{cj,cl} delta delta = 4 * N = 12.            *)
col = col /. SUNFDelta[_, _] -> SUNN;
col = Expand[col /. SUNN -> 3];
leftover = Union[Cases[col, (SUNT|SUNTF|SUNF|SUND|SUNTrace|SUNIndex|SUNFIndex|
                            SUNDelta|SUNFDelta|IndexDelta)[___], Infinity]];
Print["leftover colour objects (should be {}): ", leftover];
Print["colour-summed LC: ", LeafCount[col]];

rules = (#[[1]] :> I0 ser[#[[2]]]) & /@ num;
val = Expand[col /. rules];
val = val /. SMP["g_s"] -> 1;              (* strip the trivial g_s^4 *)
Print["after master substitution, LC: ", LeafCount[val]];
Print["non-numeric symbols left (should be {}): ",
      Select[Union@Cases[val,_Symbol,Infinity], !NumericQ[#] && # =!= ep &]];

(* Same normalisation the radiator assemblers use (amplitude_functions_1mass.wl,
   _nab.wl): the loop measure factor I0=I above turns the bare masters into the
   amplitude convention, and MSBarFac supplies the MS-bar factor.  MSBarFac carries
   a 1/I, so without it M1.M0* comes out purely imaginary at a Euclidean point and
   2Re would spuriously vanish.  M0 is tree level and gets no such factor. *)
MSBarFac = (16 Pi^2/I) (Exp[EulerGamma]/(4 Pi))^ep;
val = val MSBarFac;

ser0 = Series[val, {ep, 0, ORD}];
lau  = Table[N[SeriesCoefficient[ser0, k], 10], {k, -2, ORD}];
Print["\nQUARK one-loop x tree interference  M1.M0*  at s=-3,t=-5,u=-2,mq^2=1:"];
Do[Print["  eps^", k, " : ", Chop[lau[[k+3]], 10^-8]], {k, -2, ORD}];
Print["\n2 Re{M1.M0*}  (the physical interference):"];
Do[Print["  eps^", k, " : ", Chop[2 Re[lau[[k+3]]], 10^-8]], {k, -2, ORD}];
(* normalise by the colour/spin/pol-summed TREE from 12ferm_tree.wl:
   |M0|^2 = 712/3 - 120 eps at this point.  The ratio 2Re<M0|M1>/|M0|^2 is the
   gauge-invariant, colour-independent quantity to compare with the SCALAR
   (squark) radiator -- their difference is the magnetic remainder. *)
tree = 712/3 - 120 ep;
rel  = Series[2 Re[val /. ep -> epr] /. epr -> ep, {ep, 0, ORD}];
Print["\ntree |M0|^2 = ", tree, "   (from 12ferm_tree.wl)"];
Print["\n2Re{M1.M0*} / |M0|^2 :"];
Do[Print["  eps^", k, " : ",
   Chop[N[SeriesCoefficient[Series[(2 Re[lau[[k+3]]]) , {ep,0,0}],0],10],10^-8]], {k,-2,ORD}];
lauN = Table[N[SeriesCoefficient[Series[val/tree, {ep,0,ORD}], k],10], {k,-2,ORD}];
Print["\nfull ratio  M1.M0*/|M0|^2  (complex):"];
Do[Print["  eps^", k, " : ", Chop[lauN[[k+3]], 10^-8]], {k,-2,ORD}];
Print["\n2Re of that ratio:"];
Do[Print["  eps^", k, " : ", Chop[2 Re[lauN[[k+3]]], 10^-8]], {k,-2,ORD}];
DumpSave["/tmp/quark_interference_value.mx", lau];
Print["\nsaved /tmp/quark_interference_value.mx"];
