(* Stage B closure, RAM-frugal.  Per term: spinsum->pol->colour->substitute NUMERIC
   masters->extract eps-coeffs as NUMBERS, accumulate numbers (tiny).  Integer
   kinematics with hierarchy mq^2 << |t| << |s| (mq^2=1 fixed, hard scales large) to
   keep arithmetic integer and expose the collinear log(mq^2/hard).  Two hard scales
   give the log slope -> compare Hoeche -8(1-z), z=s/(s+u). *)
PrependTo[$Path, FileNameJoin[{$HomeDirectory,"fc93"}]];
$LoadAddOns={"FeynArts"}; $FeynCalcStartupMessages=False; Quiet[Get["FeynCalc`"]]; $FAVerbose=0;
base="/Users/user/Library/CloudStorage/Dropbox/LogMoments/LogdPhipT/splitting_1loop";
SetDirectory[base]; fcr=FileNameJoin[{base,"for_stefan_1mass","feyncalc_reproduction"}];
Get[FileNameJoin[{fcr,"colour_singlet.wl"}]];
modelM=FileNameJoin[{fcr,"SQCDBGF_massive","SQCDBGF"}];
SetMandelstam[s,t,u,-p,k2,k3,k1,Q,0,mq,mq]; ep=Global`Epsilon; ORD=0;
mu2=1; cG=(4 Pi)^ep Gamma[1+ep]Gamma[1-ep]^2/Gamma[1-2 ep];
hypE[zz_]:=Hypergeometric2F1[-ep,-ep,1-ep,N[zz,40]];
boxF[a_,b_,q_]:=Module[{cc=q-a-b},2 I cG/ep^2 (1/(a b))(
   (-(b+cc)/(a (b/mu2)))^ep hypE[1-b/(b+cc)]+(-(a+cc)/(a (b/mu2)))^ep hypE[1-a/(a+cc)]
 - (-(a+cc)(b+cc)/((a+b+cc) a (b/mu2)))^ep hypE[1-(a b)/((a+cc)(b+cc))])];
Bf[x_]:=I cG/(ep(1-2 ep))(-mu2/x)^ep; Cf1[x_]:=I cG/ep^2 (1/x)(-mu2/x)^ep;
Cf2[x_,y_]:=I cG/ep^2 (1/(x-y))((-mu2/x)^ep-(-mu2/y)^ep);
A0f[m2_]:=I cG (mu2/m2)^ep m2/(ep(1-ep)); B0mm0[m2_]:=I cG (mu2/m2)^ep/ep; B0os[m2_]:=I cG (mu2/m2)^ep/(ep(1-2 ep));
mval[sv_,tv_,uv_,m2_]:=Module[{q=sv+tv+uv-2 m2},{
  A0[m2]->A0f[m2],B0[0,m2,m2]->B0mm0[m2],B0[m2,0,m2]->B0os[m2],
  B0[sv,m2,m2]->Bf[sv],B0[tv,0,m2]->Bf[tv],B0[uv,0,m2]->Bf[uv],B0[q,m2,m2]->Bf[q],
  C0[0,m2,tv,0,0,m2]->Cf1[tv],C0[0,m2,tv,m2,m2,0]->Cf1[tv],C0[0,m2,uv,0,0,m2]->Cf1[uv],
  C0[0,m2,uv,m2,m2,0]->Cf1[uv],C0[m2,m2,sv,m2,0,m2]->Cf1[sv],
  C0[0,sv,q,m2,m2,m2]->Cf2[sv,q],C0[m2,tv,q,m2,0,m2]->Cf2[tv,q],C0[m2,uv,q,m2,0,m2]->Cf2[uv,q],
  D0[0,uv,m2,sv,m2,q,m2,m2,0,m2]->boxF[uv,sv,q],D0[m2,tv,q,sv,0,m2,m2,0,m2,m2]->boxF[sv,tv,q],
  D0[0,uv,q,tv,m2,m2,0,0,m2,m2]->boxF[uv,tv,q],D0[m2,0,m2,q,tv,uv,m2,0,0,m2]->boxF[tv,uv,q],
  D0[m2,tv,m2,uv,0,q,0,m2,0,m2]->boxF[tv,uv,q]}];
tree[out_,excl_]:=Total[FCFAConvert[CreateFeynAmp[InsertFields[CreateTopologies[0,1->3],
  {S[1,{ck,cj}]}->out,InsertionLevel->{Classes},Model->modelM,GenericModel->modelM,ExcludeParticles->excl]],
  IncomingMomenta->{p},OutgoingMomenta->{k1,k2,k3},UndoChiralSplittings->True,ChangeDimension->D,
  List->True,LoopMomenta->{},SMP->True,Contract->True,DropSumOver->True]]//FeynAmpDenominatorExplicit//ScalarProductExpand;
M0q=tree[{F[3,{1,ci}],V[50,{ca}],-F[3,{1,cl}]},{S[1],S[14],S[13],F[10]}];
M0sq=tree[{S[13,{1,ci}],V[50,{ca}],-S[13,{1,cl}]},{S[1],S[14],F[10]}];
Get["/tmp/ferm_reduced_massive.mx"]; M1q=tot;
Get["/tmp/sq_reduced_massive.mx"]; M1sq=tot;
Print["loaded. terms quark=",Length[M1q]," squark=",Length[M1sq]," memInUse=",N[MemoryInUse[]/10^9,3],"GB"];
coefs[expr_]:=Table[N[SeriesCoefficient[Series[expr,{ep,0,ORD}],k],10],{k,-2,ORD}];
frugal[M1_,M0_,ferm_,kin_,m2_,sv_,tv_,uv_]:=Module[{terms,M0c,acc,i,ti,ok},
  M0c=ComplexConjugate[M0/.kin]; terms=If[Head[M1]===Plus,List@@M1,{M1}];
  acc=ConstantArray[0,ORD+3];
  Do[ti=singletM[(terms[[i]]/.kin)] singletMbar[M0c];
     ti=If[ferm,DiracSimplify[FermionSpinSum[ti]],ti];
     ti=Contract[DoPolarizationSums[ti,k2,0]]/.D->4-2ep;
     ti=colourContract[Expand[ScalarProductExpand[ti]]]/3/.SMP["g_s"]->1;
     ti=ti/.mval[sv,tv,uv,m2];
     acc=acc+coefs[ti];
     Clear[ti]; If[Mod[i,10]==0,ClearSystemCache[];Share[]],{i,Length[terms]}];
  acc];
run[sv_,tv_,uv_,m2_]:=Module[{cq,cs},
  cq=frugal[M1q,M0q,True,{s->sv,t->tv,u->uv,mq->Sqrt[m2],mq^2->m2},m2,sv,tv,uv];
  Print["  quark done mem=",N[MemoryInUse[]/10^9,3],"GB"];
  cs=frugal[M1sq,M0sq,False,{s->sv,t->tv,u->uv,mq->Sqrt[m2],mq^2->m2},m2,sv,tv,uv];
  Print["  squark done"];
  cq-cs];
Print["target: eps^-2 coeff of collinear log slope should give -8(1-z)=-4 at z=1/2"];
Do[Module[{r=run[P[[1]],P[[2]],P[[3]],1]},
   Print[">> (s,t,u)=",P,"  z=",N[P[[1]]/(P[[1]]+P[[3]]),3],
         "  R^mag eps^-2=",Chop[r[[1]],10^-4],"  eps^-1=",Chop[r[[2]],10^-4]]],
 {P,{{-2000,-100,-2000},{-20000,-1000,-20000}}}];
