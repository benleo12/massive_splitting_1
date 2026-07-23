(* 04mm: TWO-mass version of 04m_abelian_generate_reduce_massive.wl.
   Same 6 abelian diagrams + same colour projections, but the two external
   squarks are DIFFERENT generations -> DIFFERENT masses:
     - k1 = squark      gen 1, mass m1   (S[13,{1,ci}])
     - k3 = antisquark  gen 2, mass m2   (S[13,{2,ck}])
     - k2 = gluon, massless
   Model SQCDBGF_2mass: MQU[1]=m1, MQU[2]=m2, and the S[1]-Su-Su* vertex has the
   generation IndexDelta REMOVED (flavour-off-diagonal current), so the squark
   line carries m1 on the k1 side of the S[1] vertex and m2 on the k3 side.
   Constraint becomes  s+t+u = Q^2 + m1^2 + m2^2.
   Output: ampRaw_2mass.txt + master_list_2mass.txt; sanity m2->m1 collapses to
   the 1-mass master set, m1=m2->0 to the massless set. *)
PrependTo[$Path, "/tmp/fc93"]; $LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
Print["FeynCalc ", FeynCalc`$FeynCalcVersion];
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_2mass","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,m2,m1];   (* k1=squark m1, k3=antisquark m2, k2=gluon 0 *)
CF=(SUNN^2-1)/(2 SUNN); CA=SUNN;
TaTAbaRule := {SUNTF[{a_},i_,j_](SUNTF[{b_,a_},k_,l_]+SUNTF[{a_,b_},k_,l_]) ->
   SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2)+(-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]};
proj1113 := {SUNTF[{a_},k_,m_]SUNTF[{a_},i_,j_]SUNTF[{b_},m_,l_] ->
   (1/2)((-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]+SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2))};
proj510 := {SUNTF[{a_},i_,j_]SUNTF[{b_,a_},k_,l_] ->
   (1/2)((-SUNN^(-1))SUNTF[{b},k,l]SUNFDelta[i,j]+SUNTF[{b},i,l](SUNFDelta[k,j]/2)+SUNTF[{b},k,j](SUNFDelta[i,l]/2))};

diagsl2=InsertFields[CreateTopologies[1,1->3,ExcludeTopologies->{Tadpoles}],
  {S[1,{cl,cj}]}->{S[13,{1,ci}],V[50,{ca}],S[13,{2,ck}]}, InsertionLevel->{Classes},
  Model->modelM, GenericModel->modelM, ExcludeParticles->{S[1],S[14],F[_]}];
Print["#diagrams generated: ", Length[DiagramExtract[diagsl2, All]]];
conv[idx_,lm_]:=FCFAConvert[CreateFeynAmp[DiagramExtract[diagsl2,{idx}]],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,
  ChangeDimension->D,List->False,LoopMomenta->lm,SMP->True,Contract->True,DropSumOver->True];
polk2=Pair[Momentum[k2,D],Momentum[Polarization[k2,-I],D]];

ampsl12$0   = (conv[1,{l}] /. TaTAbaRule) + ((conv[2,{l}] /. {cj->cl,cl->cj}) /. TaTAbaRule);
ampsal1113$0= Simplify[ScalarProductExpand[{
   conv[9,{-l+k1+k3}] /. proj1113, conv[11,{l+k2}] /. proj1113 }] /. {polk2->0}];
ampsal510$0 = { SUNSimplify[conv[5,{l}]] /. {cl->cj,cj->cl} /. proj510,
                SUNSimplify[conv[8,{l}]] /. proj510 };
Print["amp[0] LCs: ", LeafCount/@{ampsl12$0, Total@Flatten@{ampsal1113$0}, Total@Flatten@{ampsal510$0}}];

tid[a_]:=TID[FCE[a],l,UsePaVeBasis->True,ToPaVe->True]/.{polk2->0};
red[a_]:=PaVeReduce[FeynAmpDenominatorExplicit[a],PaVeAutoReduce->True];
Print["TID+reduce ampsl12..."];   a12  = red[tid[ampsl12$0]];
Print["TID+reduce ampsal1113..."];a1113= red[tid[Total@Flatten@{ampsal1113$0}]];
Print["TID+reduce ampsal510..."]; a510 = red[tid[Total@Flatten@{ampsal510$0}]];

sum = FCE[a12 + a1113 + a510];
polk1=SPD[k1,Polarization[k2,-I]];
polk3=SPD[k3,Polarization[k2,-I]];
abK1 = Coefficient[sum /. {polk3->0}, polk1];
Print["2-MASS abelian k1-coeff LC: ", LeafCount[abK1]];
Print["m1 present? ", !FreeQ[abK1, m1], "   m2 present? ", !FreeQ[abK1, m2]];
Print["B0 args: ", Union[Cases[abK1,B0[x__]:>{x},Infinity]]];
Print["C0 args: ", Union[Cases[abK1,C0[x__]:>{x},Infinity]]];
Print["D0 args: ", Union[Cases[abK1,D0[x__]:>{x},Infinity]]];
DumpSave["/tmp/ab_coeffK1_2mass.mx", abK1];

(* freeze portable text *)
out2m=FileNameJoin[{base,"for_stefan_2mass"}]; If[!DirectoryQ[out2m],CreateDirectory[out2m]];
Export[FileNameJoin[{out2m,"ampRaw_2mass.txt"}], ToString[abK1, InputForm], "Text"];
Print["saved ", FileNameJoin[{out2m,"ampRaw_2mass.txt"}]];

(* sanity: m2->m1 should collapse to the 1-mass master set (single mass) *)
Print["--- m2->m1 (should match 1-mass master set) ---"];
abK1m1 = abK1 /. m2->m1;
Print["B0 args (m2->m1): ", Union[Cases[abK1m1,B0[x__]:>{x},Infinity]]];
Print["C0 args (m2->m1): ", Union[Cases[abK1m1,C0[x__]:>{x},Infinity]]];
Print["D0 args (m2->m1): ", Union[Cases[abK1m1,D0[x__]:>{x},Infinity]]];
