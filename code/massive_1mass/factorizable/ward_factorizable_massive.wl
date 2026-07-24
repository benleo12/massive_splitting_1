(* TASK 1: massive Ward identity of the factorizable (Fig.2) set.
   Compute at one loop, with the squark MASS kept:
     (i)  the squark self-energy  Sigma(p^2)
     (ii) the squark-squark-gluon vertex correction Gamma^mu(q,p,k), q=p+k,
          with BOTH scalar legs off-shell (p^2, q^2 arbitrary), gluon on-shell.
   Check the naive (background-field) Ward identity
        k_mu Gamma^mu(q,p,k)  =  g_s T^a [ Sigma(q^2) - Sigma(p^2) ] ,
   which is what guarantees the abelian part of the vertex corrections cancels the
   self-energies and the factorizable set is gauge invariant with m =/= 0. *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False;
Quiet[Get["FeynCalc`"]]; $FAVerbose=0; $KeepLogDivergentScalelessIntegrals=False;
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base];
modelM=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction","SQCDBGF_massive","SQCDBGF"}];

(* ---------- (i) self-energy: S[13] -> S[13] ---------- *)
ClearAll[p,k,q]; 
FCClearScalarProducts[];
ScalarProduct[p,p]=psq;
dSE=InsertFields[CreateTopologies[1,1->1,ExcludeTopologies->{Tadpoles}],
  {S[13,{1,ci}]}->{S[13,{1,cj}]},InsertionLevel->{Classes},
  Model->modelM,GenericModel->modelM,ExcludeParticles->{S[1],S[14],F[_]}];
Print["self-energy diagrams: ",Length[CreateFeynAmp[dSE]]];
se=FCFAConvert[CreateFeynAmp[dSE],IncomingMomenta->{p},OutgoingMomenta->{p},
  UndoChiralSplittings->True,ChangeDimension->D,List->False,LoopMomenta->{l},
  SMP->True,Contract->True,DropSumOver->True];
se=TID[FCE[se],l,UsePaVeBasis->True,ToPaVe->True];
se=PaVeReduce[FeynAmpDenominatorExplicit[se],PaVeAutoReduce->True];
se=SUNSimplify[se,SUNNToCACF->True]//Simplify;
Print["Sigma(p^2) = ",InputForm[se/.SMP["g_s"]->gs]];
DumpSave["/tmp/se_massive.mx",se];

(* ---------- (ii) vertex: S[13] -> S[13] V[50], legs off-shell ---------- *)
FCClearScalarProducts[];
ScalarProduct[q,q]=qsq; ScalarProduct[p,p]=psq; ScalarProduct[k,k]=0;
ScalarProduct[p,k]=(qsq-psq)/2; ScalarProduct[q,p]=(qsq+psq)/2; ScalarProduct[q,k]=(qsq-psq)/2;
dV=InsertFields[CreateTopologies[1,1->2,ExcludeTopologies->{Tadpoles,WFCorrections}],
  {S[13,{1,ci}]}->{S[13,{1,cj}],V[50,{ca}]},InsertionLevel->{Classes},
  Model->modelM,GenericModel->modelM,ExcludeParticles->{S[1],S[14],F[_]}];
Print["vertex diagrams: ",Length[CreateFeynAmp[dV]]];
vx=FCFAConvert[CreateFeynAmp[dV],IncomingMomenta->{q},OutgoingMomenta->{p,k},
  UndoChiralSplittings->True,ChangeDimension->D,List->False,LoopMomenta->{l},
  SMP->True,Contract->True,DropSumOver->True];
(* contract with the gluon momentum k BEFORE reduction: replace the polarisation by k *)
vxk=vx/.{Momentum[Polarization[k,-I],dd_]:>Momentum[k,dd],
         Momentum[Polarization[k,-I]]:>Momentum[k]};
vxk=Contract[vxk];
vxk=TID[FCE[vxk],l,UsePaVeBasis->True,ToPaVe->True];
vxk=PaVeReduce[FeynAmpDenominatorExplicit[vxk],PaVeAutoReduce->True];
vxk=SUNSimplify[vxk,SUNNToCACF->True]//Simplify;
Print["k.Gamma = ",InputForm[vxk/.SMP["g_s"]->gs]];
DumpSave["/tmp/vxk_massive.mx",vxk];

(* ---------- Ward check: k.Gamma vs Sigma(q^2)-Sigma(p^2) ---------- *)
Get["/tmp/se_massive.mx"];
sig[x_]:=se/.psq->x;
diff=Simplify[sig[qsq]-sig[psq]];
Print["Sigma(q2)-Sigma(p2) = ",InputForm[diff/.SMP["g_s"]->gs]];
(* the WI: k.Gamma = g_s T^a (Sigma(q2)-Sigma(p2)) up to the colour factor:
   strip colour on both sides and compare ratios numerically at random points *)
stripC[x_]:=x/.{SUNTF[__]->1,SUNT[__]->1,SUNDelta[__]->1,SUNFDelta[__]->1,IndexDelta[__]->1,SUNN->3,CA->3,CF->4/3};
r1=stripC[vxk]; r2=stripC[diff]/.SMP["g_s"]->SMP["g_s"];
Do[Module[{vals={qsq->pt[[1]],psq->pt[[2]],mq->1,SMP["g_s"]->1},a,b},
  a=Simplify[r1/.vals]; b=Simplify[r2/.vals];
  Print["  at (q2,p2)=",pt,":  k.Gamma/(Sig(q2)-Sig(p2)) = ",
        Simplify[a/b]]],
 {pt,{{-3,-7},{-5,-2},{-11,-4}}}];
Print["WARD_FACT_DONE"];
