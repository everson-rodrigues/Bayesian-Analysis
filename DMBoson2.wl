(* ::Package:: *)

(*This one doesnt calculate the \[Rho]t, only use previously found values*)

BeginPackage["DMBoson2`"]
DMRoutine::usage="DMRoutine[\[Rho]cii]";
Begin["Private`"]
(*fDM_,MDM_,csd_,c\[Nu]d_*)
DMRoutine[\[Rho]cii_,\[Rho]cff_,fDMi_,MDMii_,c\[Nu]dii_,csdii_,SRi_,step\[Rho]ti_]:=Block[ {\[Rho]f,Print},
AC=8;PR=6;IT=125;
Clear[\[Omega]0,\[Rho]b0,\[Sigma]];
Clear[ \[Rho],\[Rho]0,\[Rho]n,\[Rho]p,\[Rho]s,\[Rho]e,\[Rho]2,\[Rho]b0,pp];
Clear[\[Mu]p,\[Mu]n,\[Mu]pi,\[Mu]ni,\[Mu]e];
Clear[g\[Sigma],g\[Omega],g\[Rho],m\[Sigma],m\[Omega],m\[Rho]];
Clear[\[Beta], T,\[Nu] ,\[Gamma],I\[Rho],kf,Ms,Mn,Model,\[Epsilon],P];
Clear[\[Rho]sDM,kDM,\[Epsilon]teDM,PteDM,MsDM,\[Rho]DM];
Clear[\[Phi]n,\[Phi]p,\[CapitalDelta]n,\[CapitalDelta]p];

(*Integral form of analytical expressions*)
(*CompoundExpression[]*)
\[Rho]sp[Ms_,\[Rho]pn_]=0.025330295910584444` Ms \[Gamma] (kf[\[Rho]pn] Sqrt[Ms^2+kf[\[Rho]pn]^2]+Ms^2 (-1.` Log[Ms]+Log[-1.` kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2]]));
\[Rho]sn[Ms_,\[Rho]pn_]=0.025330295910584444` Ms \[Gamma] (kf[\[Rho]pn] Sqrt[Ms^2+kf[\[Rho]pn]^2]+Ms^2 (-1.` Log[Ms]+Log[-1.` kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2]]));
\[Rho]spSR[Ms_,\[Rho]pn_,\[CapitalDelta]np_,\[Phi]np_,Cnp_]=(Cnp \[Gamma] Sqrt[1+Ms^2/kf[\[Rho]pn]^2] kf[\[Rho]pn]^4)/(2 Ms \[Pi]^2)-(Cnp \[Gamma] Sqrt[1+Ms^2/(\[Phi]np^2 kf[\[Rho]pn]^2)] kf[\[Rho]pn]^4)/(2 Ms \[Pi]^2)+(Ms \[Gamma] \[CapitalDelta]np kf[\[Rho]pn] Sqrt[Ms^2+kf[\[Rho]pn]^2])/(4 \[Pi]^2)+(Ms^3 \[Gamma] \[CapitalDelta]np Log[(-kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2])/Ms])/(4 \[Pi]^2);
\[Rho]snSR[Ms_,\[Rho]pn_,\[CapitalDelta]np_,\[Phi]np_,Cnp_]=(Cnp \[Gamma] Sqrt[1+Ms^2/kf[\[Rho]pn]^2] kf[\[Rho]pn]^4)/(2 Ms \[Pi]^2)-(Cnp \[Gamma] Sqrt[1+Ms^2/(\[Phi]np^2 kf[\[Rho]pn]^2)] kf[\[Rho]pn]^4)/(2 Ms \[Pi]^2)+(Ms \[Gamma] \[CapitalDelta]np kf[\[Rho]pn] Sqrt[Ms^2+kf[\[Rho]pn]^2])/(4 \[Pi]^2)+(Ms^3 \[Gamma] \[CapitalDelta]np Log[(-kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2])/Ms])/(4 \[Pi]^2);

(******Pressure******)
Pkin[Ms_,\[Rho]pn_]=1/(6 \[Pi]^2) \[Gamma] (-0.375` Ms^2 kf[\[Rho]pn] Sqrt[Ms^2+kf[\[Rho]pn]^2]+0.25` kf[\[Rho]pn]^3 Sqrt[Ms^2+kf[\[Rho]pn]^2]+Ms^4 (0.375` Log[Ms]-0.375` Log[-1.` kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2]]));

PkinSR[Ms_,\[Rho]pn_,\[CapitalDelta]np_,\[Phi]np_,Cnp_]=-((Ms^2 \[Gamma] \[CapitalDelta]np kf[\[Rho]pn] Sqrt[Ms^2+kf[\[Rho]pn]^2])/(16 \[Pi]^2))+(\[Gamma] \[CapitalDelta]np kf[\[Rho]pn]^3 Sqrt[Ms^2+kf[\[Rho]pn]^2])/(24 \[Pi]^2)+(Ms^4 \[Gamma] \[CapitalDelta]np Log[(kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2])/Ms])/(16 \[Pi]^2)+1/(12 \[Pi]^2) Cnp \[Gamma] kf[\[Rho]pn]^4 Log[-(((-kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2]) (\[Phi]np kf[\[Rho]pn]+Sqrt[Ms^2+\[Phi]np^2 kf[\[Rho]pn]^2]))/((kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2]) (\[Phi]np kf[\[Rho]pn]-Sqrt[Ms^2+\[Phi]np^2 kf[\[Rho]pn]^2])))];

(******Energy Density******)
\[Epsilon]kin[Ms_,\[Rho]pn_]=1/(16 \[Pi]^2) \[Gamma] (Ms^2 kf[\[Rho]pn] Sqrt[Ms^2+kf[\[Rho]pn]^2]+2 kf[\[Rho]pn]^3 Sqrt[Ms^2+kf[\[Rho]pn]^2]+Ms^4 (-Log[Ms]+Log[-kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2]]));

\[Epsilon]kinSR[Ms_,\[Rho]pn_,\[CapitalDelta]np_,\[Phi]np_,Cnp_]=(Ms^2 \[Gamma] \[CapitalDelta]np kf[\[Rho]pn] Sqrt[Ms^2+kf[\[Rho]pn]^2])/(16 \[Pi]^2)+(Cnp \[Gamma] kf[\[Rho]pn]^3 Sqrt[Ms^2+kf[\[Rho]pn]^2])/(2 \[Pi]^2)+(\[Gamma] \[CapitalDelta]np kf[\[Rho]pn]^3 Sqrt[Ms^2+kf[\[Rho]pn]^2])/(8 \[Pi]^2)-(Cnp \[Gamma] kf[\[Rho]pn]^3 Sqrt[Ms^2+\[Phi]np^2 kf[\[Rho]pn]^2])/(2 \[Pi]^2 \[Phi]np)-(Ms^4 \[Gamma] \[CapitalDelta]np Log[Ms])/(16 \[Pi]^2)+(Ms^4 \[Gamma] \[CapitalDelta]np Log[-kf[\[Rho]pn]+Sqrt[Ms^2+kf[\[Rho]pn]^2]])/(16 \[Pi]^2)-(Cnp \[Gamma] kf[\[Rho]pn]^4 Log[kf[\[Rho]pn]/Ms+Sqrt[1+kf[\[Rho]pn]^2/Ms^2]])/(2 \[Pi]^2)+(Cnp \[Gamma] kf[\[Rho]pn]^4 Log[(\[Phi]np kf[\[Rho]pn])/Ms+Sqrt[1+(\[Phi]np^2 kf[\[Rho]pn]^2)/Ms^2]])/(2 \[Pi]^2);



(*----------------------------------------Total energy and pressure--------------------------------------------*)
(*No muons*)
Ptb[\[Rho]_,\[Rho]e_,Ms_,pp_]:=P[\[Rho],Ms,pp]+\[Mu]e[\[Rho]e,Ms]^4/(12\[Pi]^2);
\[Epsilon]tb[\[Rho]_,\[Rho]e_,Ms_,pp_]:=\[Epsilon][\[Rho],Ms,pp]+\[Mu]e[\[Rho]e,Ms]^4/(4\[Pi]^2);
(*Muons*)
Pte[\[Rho]_,\[Rho]e_,Ms_,pp_]:=P[\[Rho],Ms,pp]+\[Mu]e[\[Rho]e,Ms]^4/(12\[Pi]^2)+(1/(3\[Pi]^2))NIntegrate[x^4/(x^2+m\[Mu]^2)^(1/2),{x,0,\[Sqrt](\[Mu]e[\[Rho]e,Ms]^2- m\[Mu]^2)}];
\[Epsilon]te[\[Rho]_,\[Rho]e_,Ms_,pp_]:=\[Epsilon][\[Rho],Ms,pp]+\[Mu]e[\[Rho]e,Ms]^4/(4\[Pi]^2)+(1/\[Pi]^2)NIntegrate[x^2 (x^2+m\[Mu]^2)^(1/2),{x,0,\[Sqrt](\[Mu]e[\[Rho]e,Ms]^2- m\[Mu]^2)}];

(*#############################################################################################################################*)
(*Find the transition \[Rho]*)
\[Rho]tRoutine[\[Rho]ci_,\[Rho]cf_,fDM_,c\[Nu]di_,csdi_,SR_,step\[Rho]t_,AC_,PR_,IT_]:=Block[{\[Rho]f},VtherNum[ppinp_,\[Rho]inp_]:=Block[{},Clear[\[Rho],el];i=j=1;stepd=0.00001`;TempEM={{el[1][1],el[1][2],el[1][3]},{el[2][1],el[2][2],el[2][3]},{el[3][1],el[3][2],el[3][3]}};For[pp=ppinp-stepd,pp<=ppinp+stepd,pp=pp+stepd,For[\[Rho]=\[Rho]inp-stepd,\[Rho]<=\[Rho]inp+stepd,\[Rho]=\[Rho]+stepd,eq1:=(m\[Rho]^2 \[Rho]b0-1/2 (g\[Rho] (2 pp-1)) \[Rho]+(\[Alpha]3L ((g\[Omega] g\[Rho]) \[Omega]0)^2) \[Rho]b0)/Mn^3;eq2:=(m\[Omega]^2 \[Omega]0-g\[Omega] \[Rho]+(Cc g\[Omega]) (g\[Omega] \[Omega]0)^3+(\[Alpha]3L ((g\[Omega] g\[Rho]) \[Rho]b0)^2) \[Omega]0)/Mn^3;eq3:=((m\[Sigma]^2 (Mn-Ms))/g\[Sigma]-g\[Sigma] \[Rho]s[\[Rho],Ms,pp]+A ((Mn-Ms)/g\[Sigma])^2+B ((Mn-Ms)/g\[Sigma])^3)/Mn^3;Ans\[Omega]\[Rho]=FindRoot[{(m\[Rho]^2 \[Rho]b0-1/2 (g\[Rho] (2 pp-1)) \[Rho]+(\[Alpha]3L ((g\[Omega] g\[Rho]) \[Omega]0)^2) \[Rho]b0)/Mn^3,(m\[Omega]^2 \[Omega]0-g\[Omega] \[Rho]+(Cc g\[Omega]) (g\[Omega] \[Omega]0)^3+(\[Alpha]3L ((g\[Omega] g\[Rho]) \[Rho]b0)^2) \[Omega]0)/Mn^3,((m\[Sigma]^2 (Mn-Ms))/g\[Sigma]-g\[Sigma] \[Rho]s[\[Rho],Ms,pp]+A ((Mn-Ms)/g\[Sigma])^2+B ((Mn-Ms)/g\[Sigma])^3)/Mn^3},{{Ms,0.001`,Mn},{\[Omega]0,0.0001`,0.1`},{\[Rho]b0,-0.05`,0.05`}},Evaluated->False,MaxIterations->IT,AccuracyGoal->AC,PrecisionGoal->PR];\[Rho]e=\[Rho]e/. Ans\[Omega]\[Rho];\[Rho]b0=\[Rho]b0/. Ans\[Omega]\[Rho];\[Omega]0=\[Omega]0/. Ans\[Omega]\[Rho];Ms=Ms/. Ans\[Omega]\[Rho];pp=pp/. Ans\[Omega]\[Rho];\[Sigma]=(Mn-Ms)/g\[Sigma];el[i][j]=\[Epsilon][\[Rho],Ms,pp]/\[Rho];i=i+1;];i=1;j=j+1;];D\[Epsilon]D\[Rho]=(TempEM[[3]][[2]]-TempEM[[1]][[2]])/(2 stepd);D\[Epsilon]2D\[Rho]2=(TempEM[[1]][[2]]-2 TempEM[[2]][[2]]+TempEM[[3]][[2]])/stepd^2;D\[Epsilon]Dy=(TempEM[[2]][[3]]-TempEM[[2]][[1]])/(2 stepd);D\[Epsilon]2Dy2=(TempEM[[2]][[1]]-2 TempEM[[2]][[2]]+TempEM[[2]][[3]])/stepd^2;D\[Epsilon]2DyD\[Rho]=(TempEM[[1]][[1]]-TempEM[[1]][[3]]-TempEM[[3]][[1]]+TempEM[[3]][[3]])/(4 stepd^2);Vthermn=(2 \[Rho]) D\[Epsilon]D\[Rho]+\[Rho]^2 D\[Epsilon]2D\[Rho]2-(\[Rho] D\[Epsilon]2DyD\[Rho])^2/D\[Epsilon]2Dy2;AppendTo[D\[Epsilon]D\[Rho]List,D\[Epsilon]D\[Rho]];AppendTo[D\[Epsilon]2D\[Rho]2List,D\[Epsilon]2D\[Rho]2];AppendTo[D\[Epsilon]DyList,D\[Epsilon]Dy];AppendTo[D\[Epsilon]2Dy2List,D\[Epsilon]2Dy2];AppendTo[D\[Epsilon]2DyD\[Rho]List,D\[Epsilon]2DyD\[Rho]];];D\[Epsilon]D\[Rho]List={};D\[Epsilon]2D\[Rho]2List={};D\[Epsilon]DyList={};D\[Epsilon]2Dy2List={};D\[Epsilon]2DyD\[Rho]List={};If[SR==1,\[Rho]ti=0.3` n0;\[Rho]tf=0.4` n0,\[Rho]ti=0.5` n0;\[Rho]tf=0.6` n0];For[\[Rho]t=\[Rho]tf,\[Rho]t>=\[Rho]ti,\[Rho]t=\[Rho]t-step\[Rho]t,eq1:=(m\[Rho]^2 \[Rho]b0-1/2 (g\[Rho] (2 pp-1)) \[Rho]t+(\[Alpha]3L ((g\[Omega] g\[Rho]) \[Omega]0)^2) \[Rho]b0)/Mn^3;eq2:=(m\[Omega]^2 \[Omega]0-g\[Omega] \[Rho]t+(Cc g\[Omega]) (g\[Omega] \[Omega]0)^3+(\[Alpha]3L ((g\[Omega] g\[Rho]) \[Rho]b0)^2) \[Omega]0)/Mn^3;eq3:=\[Mu]n[\[Rho]t,Ms,pp]/Mn-\[Mu]p[\[Rho]t,Ms,pp]/Mn-\[Mu]e[\[Rho]e,Ms]/Mn;eq4:=(pp \[Rho]t)/n0-\[Rho]e/n0;eq5:=((m\[Sigma]^2 (Mn-Ms))/g\[Sigma]-g\[Sigma] \[Rho]s[\[Rho]t,Ms,pp]+A ((Mn-Ms)/g\[Sigma])^2+B ((Mn-Ms)/g\[Sigma])^3)/Mn^3;\[Sigma]=(Mn-Ms)/g\[Sigma];Msc=0.5` Mn;Msle=0.15` Mn;Msld=Mse;ppc=ppe;pple=0.0003`;ppld=ppe;\[Rho]ec=\[Rho]ee;\[Rho]ele=0.0001` n0;\[Rho]eld=\[Rho]ee;\[Omega]0c=\[Omega]0e;\[Omega]0le=0.0005`;\[Omega]0ld=2 \[Omega]0le;\[Rho]b0c=\[Rho]b0e;\[Rho]b0le=-0.05`;\[Rho]b0ld=0.05`;Ans\[Omega]\[Rho]t=FindRoot[{eq1,eq2,eq3,eq4,eq5},{{Ms,Msle,Msld},{pp,pple,ppld},{\[Rho]e,\[Rho]ele,\[Rho]eld},{\[Omega]0,\[Omega]0le,\[Omega]0ld},{\[Rho]b0,\[Rho]b0le,\[Rho]b0ld}},Evaluated->False,MaxIterations->IT,AccuracyGoal->AC,PrecisionGoal->PR];\[Rho]e=Re[\[Rho]e/. Ans\[Omega]\[Rho]t];\[Rho]b0=Re[\[Rho]b0/. Ans\[Omega]\[Rho]t];\[Omega]0=Re[\[Omega]0/. Ans\[Omega]\[Rho]t];Ms=Re[Ms/. Ans\[Omega]\[Rho]t];pp=Re[pp/. Ans\[Omega]\[Rho]t];If[Length[VtherList]>2&&VtherList[[-2]]<0,Print["1st Equation= ",eq1];Print["2nd Equation= ",eq2];Print["3rd Equation= ",eq3];Print["4th Equation= ",eq4];Print["5th Equation= ",eq5];];\[Sigma]=(Mn-Ms)/g\[Sigma];Block[{Print},\[Mu]numroutin[\[Rho]t,Ms,pp]];If[\[Mu]nN/\[Mu]n[\[Rho]t,Ms,pp]<0.96`||\[Mu]pN/\[Mu]p[\[Rho]t,Ms,pp]<0.96`,Break[];Print["Chemichal Potential Inconsistency"]];VtherNum[pp,\[Rho]t];If[Length[VtherList]>2&&VtherList[[-2]]<0,Print["==================================== \[Rho]=",\[Rho]t/n0," ==============================================="];Print["--------------------Chemichal Potential---------------------------"];Print[" N An= ",h \[Mu]n[\[Rho]t,Ms,pp]," |  N Num= ",h \[Mu]nN," |  Fr= ",\[Mu]nN/\[Mu]n[\[Rho]t,Ms,pp],"  |  P An= ",\[Mu]p[\[Rho]t,Ms,pp] h," | P Num= ",\[Mu]pN h," |   Fr= ",\[Mu]pN/\[Mu]p[\[Rho]t,Ms,pp]];Print["Total An= ",h ((1-pp) \[Mu]n[\[Rho]t,Ms,pp]+\[Mu]p[\[Rho]t,Ms,pp] pp)," |  Total Num= ",h ((1-pp) \[Mu]nN+\[Mu]pN pp),"  |  Euler= ",((P[\[Rho]t,Ms,pp]+\[Epsilon][\[Rho]t,Ms,pp]) h)/\[Rho]t];Print["-------------------- Vther ---------------------------"];Print["Vther= ",VtherList[[-1]]/Mn];Print["\[Rho]t= ",\[Rho]t/n0,"n0= ",\[Rho]t," \!\(\*SuperscriptBox[\(fm\), \(-3\)]\), yt= ",pp,", Pt= ",Ptt=P[\[Rho]t,Ms,pp]+\[Mu]e[\[Rho]e,Ms]^4/(12 \[Pi]^2)," \!\(\*SuperscriptBox[\(fm\), \(-4\)]\), \[Epsilon]t= ",\[Epsilon]tt=\[Epsilon][\[Rho]t,Ms,pp]+\[Mu]e[\[Rho]e,Ms]^4/(4 \[Pi]^2)," \!\(\*SuperscriptBox[\(fm\), \(-4\)]\)"];Print["Mst= ",Ms,", y= ",pp,", \[Rho]e= ",\[Rho]e,", \[Omega]0= ",\[Omega]0,", \[Rho]b0= ",\[Rho]b0];Break[]];AppendTo[VtherList,Vthermn];AppendTo[\[Rho]tList,\[Rho]t];Clear[pp,Ms,\[Rho]b0,\[Omega]0,\[Rho]e];];If[Length[Position[VtherList,-Min[Abs[VtherList]]]]==1,\[Rho]t=\[Rho]tList[[Position[VtherList,-Min[Abs[VtherList]]][[1,1]]]],\[Rho]t=\[Rho]tList[[Position[VtherList,Min[Abs[VtherList]]][[1,1]]]]];\[Rho]b0=Re[\[Rho]b0/. Ans\[Omega]\[Rho]t];\[Omega]0=Re[\[Omega]0/. Ans\[Omega]\[Rho]t];Ms=Re[Ms/. Ans\[Omega]\[Rho]t];pp=Re[pp/. Ans\[Omega]\[Rho]t];\[Rho]et=Re[\[Rho]e/. Ans\[Omega]\[Rho]t];\[Rho]b0t=Re[\[Rho]b0/. Ans\[Omega]\[Rho]t];\[Omega]0t=Re[\[Omega]0/. Ans\[Omega]\[Rho]t];Mst=Re[Ms/. Ans\[Omega]\[Rho]t];ppt=Re[pp/. Ans\[Omega]\[Rho]t];];

(*#############################################################################################################################*)
(*Zeros Finder*)
contador=0;
ZerosFinder[\[Rho]t_]:=Block[{},eq1:=(m\[Rho]^2 \[Rho]b0-1/2 (g\[Rho] (2 pp-1)) \[Rho]t+(\[Alpha]3L ((g\[Omega] g\[Rho]) \[Omega]0)^2) \[Rho]b0)/Mn^3;eq2:=(m\[Omega]^2 \[Omega]0-g\[Omega] \[Rho]t+(Cc g\[Omega]) (g\[Omega] \[Omega]0)^3+(\[Alpha]3L ((g\[Omega] g\[Rho]) \[Rho]b0)^2) \[Omega]0)/Mn^3;eq3:=\[Mu]n[\[Rho]t,Ms,pp]/Mn-\[Mu]p[\[Rho]t,Ms,pp]/Mn-\[Mu]e[\[Rho]e,Ms]/Mn;If[\[Rho]t<=(0 0.8`) n0,eq4:=(pp \[Rho]t)/n0-\[Rho]e/n0,eq4:=(pp \[Rho]t)/n0-\[Rho]e/n0-(\[Mu]e[\[Rho]e,Ms]^2-m\[Mu]^2)^(3/2)/((3 \[Pi]^2) n0)];eq5:=((m\[Sigma]^2 (Mn-Ms))/g\[Sigma]-g\[Sigma] \[Rho]s[\[Rho]t,Ms,pp]+A ((Mn-Ms)/g\[Sigma])^2+B ((Mn-Ms)/g\[Sigma])^3)/Mn^3;\[Sigma]=(Mn-Ms)/g\[Sigma];If[\[Rho]t<=0.8` n0,Msc=0.5` Mn;Msle=0.5`;Msld=Mse;ppc=ppe;pple=0.0003`;ppld=ppe;\[Rho]ec=\[Rho]ee;\[Rho]ele=0.0001` n0;\[Rho]eld=\[Rho]ee;\[Omega]0c=\[Omega]0e;\[Omega]0le=0.0005`;\[Omega]0ld=2 \[Omega]0le;\[Rho]b0c=\[Rho]b0e;\[Rho]b0le=-0.05`;\[Rho]b0ld=0.05`;,Msc=0.5` Mn;Msle=1;Msld=Mn;ppc=ppe;pple=0.09396075448574265`;ppld=0.5`;\[Rho]ec=\[Rho]ee;\[Rho]ele=0.014131273805613716`;\[Rho]eld=2 0.014131273805613716`;\[Omega]0c=\[Omega]0e;\[Omega]0le=0.10533360447998043`;\[Omega]0ld=2 0.10533360447998043`;\[Rho]b0c=\[Rho]b0e;\[Rho]b0le=-0.05`;\[Rho]b0ld=0.05`;];Ans\[Omega]\[Rho]=FindRoot[{eq1,eq2,eq3,eq4,eq5},{{Ms,Msle,Msld},{pp,pple,ppld},{\[Rho]e,\[Rho]ele,\[Rho]eld},{\[Omega]0,\[Omega]0le,\[Omega]0ld},{\[Rho]b0,\[Rho]b0le,\[Rho]b0ld}},AccuracyGoal->8,PrecisionGoal->6,Evaluated->False,MaxIterations->250];\[Rho]e=Re[\[Rho]e/. Ans\[Omega]\[Rho]];\[Rho]b0=Re[\[Rho]b0/. Ans\[Omega]\[Rho]];\[Omega]0=Re[\[Omega]0/. Ans\[Omega]\[Rho]];Ms=Re[Ms/. Ans\[Omega]\[Rho]];pp=Re[pp/. Ans\[Omega]\[Rho]];\[Sigma]=(Mn-Ms)/g\[Sigma];Print["1st Equation= ",eq1];Print["2nd Equation= ",eq2];Print["3rd Equation= ",eq3];Print["4th Equation= ",eq4];Print["5th Equation= ",eq5];If[Abs[eq1]>=1/10^6||Abs[eq1]>=1/10^6||Abs[eq2]>=1/10^6||Abs[eq3]>=1/10^6||Abs[eq4]>=1/10^6||Abs[eq5]>=1/10^6,Print[Style["############################################",Red]];Print[Style["Bad Solution at NM",Red]];notsolved=notsolved+1;]];
(*#############################################################################################################################*)
 (*\[Mu] calculator*)
\[Mu]numroutin[\[Rho]t_,Ms_,pp_]:=Block[{},\[Sigma]=(Mn-Ms)/g\[Sigma];step\[Rho]=0.01` n0;stepp=0.01` pp;dedr=(\[Epsilon][\[Rho]t+step\[Rho],Ms,pp]-\[Epsilon][\[Rho]t-step\[Rho],Ms,pp])/(2 step\[Rho]);dedy=(\[Epsilon][\[Rho]t,Ms,pp+stepp]-\[Epsilon][\[Rho]t,Ms,pp-stepp])/(2 stepp);\[Mu]pN=dedr+(((1-pp) \[Rho]t) dedy)/\[Rho]t^2;\[Mu]nN=dedr-((pp \[Rho]t) dedy)/\[Rho]t^2;Print["==================================== \[Rho]=",\[Rho]t/n0," ==============================================="];Print["Chemichal Potential"];Print[" N An= ",h \[Mu]n[\[Rho]t,Ms,pp]," |  N Num= ",h \[Mu]nN," |  Fr= ",\[Mu]nN/\[Mu]n[\[Rho]t,Ms,pp],"  |  P An= ",\[Mu]p[\[Rho]t,Ms,pp] h," | P Num= ",\[Mu]pN h," |   Fr= ",\[Mu]pN/\[Mu]p[\[Rho]t,Ms,pp]];Print["Total An= ",h ((1-pp) \[Mu]n[\[Rho]t,Ms,pp]+\[Mu]p[\[Rho]t,Ms,pp] pp)," |  Total Num= ",h ((1-pp) \[Mu]nN+\[Mu]pN pp),"  |  Euler= ",((P[\[Rho]t,Ms,pp]+\[Epsilon][\[Rho]t,Ms,pp]) h)/\[Rho]t];Print["================================================================================================================="];];

(*#############################################################################################################################*)
 (*Zeros DM calculator*)
ZerosDMFinder[\[Epsilon]tDMc_,csdii2_]:=Block[{},
\[Rho]DM=( -MDM+Sqrt[MDM^2+2 csdii2^2 \[Epsilon]tDMc]) /csdii2^2;];

(*#############################################################################################################################*)(*#############################################################################################################################*)(*#############################################################################################################################*)
(*#############################################################################################################################*)(*#############################################################################################################################*)(*#############################################################################################################################*)


contador=0;
(*SRC Parameters*)
\[Phi]n[pp_]:=\[Phi]0(1+\[Phi]1(1-2pp));
\[Phi]p[pp_]:=\[Phi]0(1-\[Phi]1(1-2pp));

Cn[pp_]:=C0(1+C1(1-2pp));
Cp[pp_]:=C0(1-C1(1-2pp));

\[CapitalDelta]p[pp_]:=(1-3Cp[pp](1-1/\[Phi]p[pp]));
\[CapitalDelta]n[pp_]:=(1-3Cn[pp](1-1/\[Phi]n[pp]));

If[SRi==0,
\[Phi]0=1;\[Phi]1=0;C0=0.161;C1=-0.25;
(*Pressure and energy density*)
\[Epsilon][\[Rho]_,Ms_,pp_]:=(1/2)(m\[Sigma] \[Sigma])^2 +(A/3)\[Sigma]^3+(B/4)\[Sigma]^4-(1/2)(m\[Omega] \[Omega]0)^2-(Cc/4)(g\[Omega] \[Omega]0)^4-(1/2)(m\[Rho] \[Rho]b0)^2+g\[Omega] \[Omega]0 (pp \[Rho]+(1-pp)\[Rho](*=\[Rho]*))+(1/2)(g\[Rho] \[Rho]b0 (pp \[Rho]-(1-pp)\[Rho])(*=\[Rho]3*)) - g\[Sigma] g\[Omega]^2 \[Sigma] \[Omega]0^2 (\[Alpha]1+1/2\[Alpha]1L g\[Sigma] \[Sigma])- g\[Sigma] g\[Rho]^2 \[Sigma] \[Rho]b0^2 (\[Alpha]2+1/2\[Alpha]2L g\[Sigma] \[Sigma])-(1/2)\[Alpha]3L (g\[Omega] g\[Rho] \[Omega]0 \[Rho]b0)^2+\[Epsilon]kinSR[Ms,pp \[Rho],\[CapitalDelta]p[pp],\[Phi]p[pp],Cp[pp]]+\[Epsilon]kinSR[Ms,(1-pp)\[Rho],\[CapitalDelta]n[pp],\[Phi]n[pp],Cn[pp]];

P[\[Rho]_,Ms_,pp_]:=-(1/2) (m\[Sigma] \[Sigma])^2 -(A/3)\[Sigma]^3-(B/4)\[Sigma]^4 +(1/2)(m\[Omega] \[Omega]0)^2+(Cc/4)(g\[Omega] \[Omega]0)^4+(1/2)(m\[Rho] \[Rho]b0)^2+(1/2)\[Alpha]3L (g\[Omega] g\[Rho] \[Omega]0 \[Rho]b0)^2+ g\[Sigma] g\[Omega]^2 \[Sigma] \[Omega]0^2 (\[Alpha]1+1/2\[Alpha]1L g\[Sigma] \[Sigma])+ g\[Sigma] g\[Rho]^2 \[Sigma] \[Rho]b0^2 (\[Alpha]2+1/2\[Alpha]2L g\[Sigma] \[Sigma])+PkinSR[Ms,pp \[Rho],\[CapitalDelta]p[pp],\[Phi]p[pp],Cp[pp]]+PkinSR[Ms,(1-pp)\[Rho],\[CapitalDelta]n[pp],\[Phi]n[pp],Cn[pp]];

\[Rho]s[\[Rho]_,Ms_,pp_]:=\[Rho]spSR[Ms,pp \[Rho],\[CapitalDelta]p[pp],\[Phi]p[pp],Cp[pp]]+\[Rho]snSR[Ms,(1-pp)\[Rho],\[CapitalDelta]n[pp],\[Phi]n[pp],Cn[pp]];
\[Mu]kinsrc[Pattern[Ms, Blank[]], Pattern[kf, Blank[]], Pattern[\[CapitalDelta], Blank[]], Pattern[\[Phi], Blank[]], Pattern[c, Blank[]]] := (3 c) (Sqrt[kf^2 + Ms^2] - Sqrt[Ms^2 + kf^2 \[Phi]^2]/\[Phi]) + ((4 c) kf) Log[(kf \[Phi] + Sqrt[Ms^2 + kf^2 \[Phi]^2])/(kf + Sqrt[kf^2 + Ms^2])] ;\[Mu]Cp[Pattern[\[Rho], Blank[]], Pattern[Ms, Blank[]], Pattern[pp, Blank[]]] := (\[Phi]0 \[Phi]1) (Cp[pp] (Sqrt[kf[pp \[Rho]]^2 + Ms^2]/\[Phi]p[pp]^2) - Cn[pp] (Sqrt[kf[\[Rho] (1 - pp)]^2 + Ms^2]/\[Phi]n[pp]^2)) + (C0 C1) (kf[pp \[Rho]]^4 (-ArcSinh[kf[pp \[Rho]]/Ms] + ArcSinh[\[Phi]p[pp] (kf[pp \[Rho]]/Ms)] + Sqrt[Ms^2 + kf[pp \[Rho]]^2]/kf[pp \[Rho]] - Sqrt[Ms^2 + \[Phi]p[pp]^2 kf[pp \[Rho]]^2]/(kf[pp \[Rho]] \[Phi]p[pp]))) - (C0 C1) (kf[(1 - pp) \[Rho]]^4 (-ArcSinh[kf[(1 - pp) \[Rho]]/Ms] + ArcSinh[\[Phi]n[pp] (kf[(1 - pp) \[Rho]]/Ms)] + Sqrt[Ms^2 + kf[(1 - pp) \[Rho]]^2]/kf[(1 - pp) \[Rho]] - Sqrt[Ms^2 + \[Phi]n[pp]^2 kf[(1 - pp) \[Rho]]^2]/(kf[(1 - pp) \[Rho]] \[Phi]n[pp]))) - ((1/8) (3 ((C0 C1) (1 - 1/\[Phi]p[pp]) + (\[Phi]0 \[Phi]1) (Cp[pp]/\[Phi]p[pp]^2)))) ((Ms^2 kf[pp \[Rho]]) Sqrt[Ms^2 + kf[pp \[Rho]]^2] + (2 kf[pp \[Rho]]^3) Sqrt[Ms^2 + kf[pp \[Rho]]^2] + Ms^4 (-Log[Ms] + Log[-kf[pp \[Rho]] + Sqrt[Ms^2 + kf[pp \[Rho]]^2]])) + ((1/8) (3 ((C0 C1) (1 - 1/\[Phi]n[pp]) + (\[Phi]0 \[Phi]1) (Cn[pp]/\[Phi]n[pp]^2)))) ((Ms^2 kf[(1 - pp) \[Rho]]) Sqrt[Ms^2 + kf[(1 - pp) \[Rho]]^2] + (2 kf[(1 - pp) \[Rho]]^3) Sqrt[Ms^2 + kf[(1 - pp) \[Rho]]^2] + Ms^4 (-Log[Ms] + Log[-kf[(1 - pp) \[Rho]] + Sqrt[Ms^2 + kf[(1 - pp) \[Rho]]^2]]));\[Mu]Cn[Pattern[\[Rho], Blank[]], Pattern[Ms, Blank[]], Pattern[pp, Blank[]]] := ((-\[Phi]0) \[Phi]1) ((Cp[pp]/\[Phi]p[pp]^2) Sqrt[kf[pp \[Rho]]^2 + Ms^2] - (Cn[pp]/\[Phi]n[pp]^2) Sqrt[kf[\[Rho] (1 - pp)]^2 + Ms^2]) - (C0 C1) (kf[pp \[Rho]]^4 (-ArcSinh[kf[pp \[Rho]]/Ms] + ArcSinh[\[Phi]p[pp] (kf[pp \[Rho]]/Ms)] + Sqrt[Ms^2 + kf[pp \[Rho]]^2]/kf[pp \[Rho]] - Sqrt[Ms^2 + \[Phi]p[pp]^2 kf[pp \[Rho]]^2]/(kf[pp \[Rho]] \[Phi]p[pp]))) + (C0 C1) (kf[(1 - pp) \[Rho]]^4 (-ArcSinh[kf[(1 - pp) \[Rho]]/Ms] + ArcSinh[\[Phi]n[pp] (kf[(1 - pp) \[Rho]]/Ms)] + Sqrt[Ms^2 + kf[(1 - pp) \[Rho]]^2]/kf[(1 - pp) \[Rho]] - Sqrt[Ms^2 + \[Phi]n[pp]^2 kf[(1 - pp) \[Rho]]^2]/(kf[(1 - pp) \[Rho]] \[Phi]n[pp]))) + (3 ((C0 C1) (1 - 1/\[Phi]p[pp]) + (\[Phi]0 \[Phi]1) (Cp[pp]/\[Phi]p[pp]^2))) ((1/8) ((Ms^2 kf[pp \[Rho]]) Sqrt[Ms^2 + kf[pp \[Rho]]^2] + (2 kf[pp \[Rho]]^3) Sqrt[Ms^2 + kf[pp \[Rho]]^2] + Ms^4 (-Log[Ms] + Log[-kf[pp \[Rho]] + Sqrt[Ms^2 + kf[pp \[Rho]]^2]]))) - (3 ((C0 C1) (1 - 1/\[Phi]n[pp]) + (\[Phi]0 \[Phi]1) (Cn[pp]/\[Phi]n[pp]^2))) ((1/8) ((Ms^2 kf[(1 - pp) \[Rho]]) Sqrt[Ms^2 + kf[(1 - pp) \[Rho]]^2] + (2 kf[(1 - pp) \[Rho]]^3) Sqrt[Ms^2 + kf[(1 - pp) \[Rho]]^2] + Ms^4 (-Log[Ms] + Log[-kf[(1 - pp) \[Rho]] + Sqrt[Ms^2 + kf[(1 - pp) \[Rho]]^2]])));\[Mu]po[Pattern[\[Rho], Blank[]], Pattern[Ms, Blank[]], Pattern[pp, Blank[]]] := Sqrt[kf[pp \[Rho]]^2 + Ms^2] \[CapitalDelta]p[pp] + \[Mu]kinsrc[Ms, kf[pp \[Rho]], \[CapitalDelta]p[pp], \[Phi]p[pp], Cp[pp]] + g\[Omega] \[Omega]0 + ((1/2) g\[Rho]) \[Rho]b0 + ((2/Pi^2) (pp (\[Rho]/\[Rho]^2))) \[Mu]Cp[\[Rho], Ms, pp];\[Mu]no[Pattern[\[Rho], Blank[]], Pattern[Ms, Blank[]], Pattern[pp, Blank[]]] := Sqrt[kf[(1 - pp) \[Rho]]^2 + Ms^2] \[CapitalDelta]n[pp] + \[Mu]kinsrc[Ms, kf[(1 - pp) \[Rho]], \[CapitalDelta]n[pp], \[Phi]n[pp], Cn[pp]] + g\[Omega] \[Omega]0 - ((1/2) g\[Rho]) \[Rho]b0 + ((2/Pi^2) ((1 - pp) (\[Rho]/\[Rho]^2))) \[Mu]Cn[\[Rho], Ms, pp];
,
(*Pressure and energy density SRC*)
\[Phi]0=2.38;\[Phi]1=-0.56;C0=0.161;C1=-0.25;

\[Epsilon][\[Rho]_,Ms_,pp_]:=(1/2)(m\[Sigma] \[Sigma])^2 +(A/3)\[Sigma]^3+(B/4)\[Sigma]^4-(1/2)(m\[Omega] \[Omega]0)^2-(Cc/4)(g\[Omega] \[Omega]0)^4-(1/2)(m\[Rho] \[Rho]b0)^2+g\[Omega] \[Omega]0 (pp \[Rho]+(1-pp)\[Rho](*=\[Rho]*))+(1/2)(g\[Rho] \[Rho]b0 (pp \[Rho]-(1-pp)\[Rho])(*=\[Rho]3*)) - g\[Sigma] g\[Omega]^2 \[Sigma] \[Omega]0^2 (\[Alpha]1+1/2\[Alpha]1L g\[Sigma] \[Sigma])- g\[Sigma] g\[Rho]^2 \[Sigma] \[Rho]b0^2 (\[Alpha]2+1/2\[Alpha]2L g\[Sigma] \[Sigma])-(1/2)\[Alpha]3L (g\[Omega] g\[Rho] \[Omega]0 \[Rho]b0)^2+\[Epsilon]kinSR[Ms,pp \[Rho],\[CapitalDelta]p[pp],\[Phi]p[pp],Cp[pp]]+\[Epsilon]kinSR[Ms,(1-pp)\[Rho],\[CapitalDelta]n[pp],\[Phi]n[pp],Cn[pp]];

P[\[Rho]_,Ms_,pp_]:=-(1/2) (m\[Sigma] \[Sigma])^2 -(A/3)\[Sigma]^3-(B/4)\[Sigma]^4 +(1/2)(m\[Omega] \[Omega]0)^2+(Cc/4)(g\[Omega] \[Omega]0)^4+(1/2)(m\[Rho] \[Rho]b0)^2+(1/2)\[Alpha]3L (g\[Omega] g\[Rho] \[Omega]0 \[Rho]b0)^2+ g\[Sigma] g\[Omega]^2 \[Sigma] \[Omega]0^2 (\[Alpha]1+1/2\[Alpha]1L g\[Sigma] \[Sigma])+ g\[Sigma] g\[Rho]^2 \[Sigma] \[Rho]b0^2 (\[Alpha]2+1/2\[Alpha]2L g\[Sigma] \[Sigma])+PkinSR[Ms,pp \[Rho],\[CapitalDelta]p[pp],\[Phi]p[pp],Cp[pp]]+PkinSR[Ms,(1-pp)\[Rho],\[CapitalDelta]n[pp],\[Phi]n[pp],Cn[pp]];

\[Rho]s[\[Rho]_,Ms_,pp_]:=\[Rho]spSR[Ms,pp \[Rho],\[CapitalDelta]p[pp],\[Phi]p[pp],Cp[pp]]+\[Rho]snSR[Ms,(1-pp)\[Rho],\[CapitalDelta]n[pp],\[Phi]n[pp],Cn[pp]];

(*############################################################################################################################*)










(*###

Part of the code omitted for confidentiality reasons.

###*)








(*############################################################################################################################*)
(*Store Normal Matter quantities*)
AppendTo[PtList[c\[Nu]dii,csdii,fDMi],Re[Ptb[\[Rho],\[Rho]e,Ms,pp]]];
AppendTo[\[Epsilon]tList[c\[Nu]dii,csdii,fDMi],Re[\[Epsilon]tb[\[Rho],\[Rho]e,Ms,pp]]];
AppendTo[\[Rho]List[c\[Nu]dii,csdii,fDMi],\[Rho]];
AppendTo[\[Mu]pList[c\[Nu]dii,csdii,fDMi],\[Mu]p[\[Rho],Ms,pp]];
AppendTo[\[Mu]nList[c\[Nu]dii,csdii,fDMi],\[Mu]n[\[Rho],Ms,pp]];
AppendTo[\[Mu]pNList[c\[Nu]dii,csdii,fDMi],\[Mu]pN];
AppendTo[\[Mu]nNList[c\[Nu]dii,csdii,fDMi],\[Mu]nN];
(*Store Dark Matter quantities*)
AppendTo[\[Rho]DMList[c\[Nu]dii,csdii,fDMi],\[Rho]DM];
AppendTo[kDMList[c\[Nu]dii,csdii,fDMi], kDM];
AppendTo[MsDMList[c\[Nu]dii,csdii,fDMi], MDM];
AppendTo[PtDMList[c\[Nu]dii,csdii,fDMi],PteDM[\[Rho]DM,MDM,kDM]];
AppendTo[\[Epsilon]tDMList[c\[Nu]dii,csdii,fDMi],\[Epsilon]teDM[\[Rho]DM,MDM,kDM]];];
,
If[\[Rho]<2n0,step\[Rho]=0.001n0];
If[\[Rho]>2n0,step\[Rho]=0.01n0];
Print[\[Rho]];
total=total+1;
ZerosFinder[\[Rho]];
\[Mu]numroutin[\[Rho],Ms,pp];
If[Re[Pte[\[Rho],\[Rho]e,Ms,pp]]>0&&Re[\[Epsilon]te[\[Rho],\[Rho]e,Ms,pp]]>0&&Im[Pte[\[Rho],\[Rho]e,Ms,pp]]<10^-8, 
(*Calculates the system of equation for DARK MATTER*)
\[Epsilon]tDM=fDMi \[Epsilon]te[\[Rho],\[Rho]e,Ms,pp] ;
ZerosDMFinder[\[Epsilon]tDM,csdii];
If[Re[PteDM[\[Rho]DM,MDM,kDM]]>0&&Re[\[Epsilon]teDM[\[Rho]DM,MDM,kDM]]>0&&Abs[Im[PteDM[\[Rho]DM,MDM,kDM]]]<10^-8&&Abs[Im[\[Epsilon]teDM[\[Rho]DM,MDM,kDM]]]<10^-8, 
(*############################################################################################################################*)
(*Store Normal Matter quantities*)
AppendTo[PtList[c\[Nu]dii,csdii,fDMi],Re[Pte[\[Rho],\[Rho]e,Ms,pp]]];
AppendTo[\[Epsilon]tList[c\[Nu]dii,csdii,fDMi],Re[\[Epsilon]te[\[Rho],\[Rho]e,Ms,pp]]];
AppendTo[\[Rho]List[c\[Nu]dii,csdii,fDMi],\[Rho]];
AppendTo[\[Mu]pList[c\[Nu]dii,csdii,fDMi],\[Mu]p[\[Rho],Ms,pp]];
AppendTo[\[Mu]nList[c\[Nu]dii,csdii,fDMi],\[Mu]n[\[Rho],Ms,pp]];
AppendTo[\[Mu]pNList[c\[Nu]dii,csdii,fDMi],\[Mu]pN];
AppendTo[\[Mu]nNList[c\[Nu]dii,csdii,fDMi],\[Mu]nN];
(*Store Dark Matter quantities*)
AppendTo[\[Rho]DMList[c\[Nu]dii,csdii,fDMi],\[Rho]DM];
AppendTo[kDMList[c\[Nu]dii,csdii,fDMi], kDM];
AppendTo[MsDMList[c\[Nu]dii,csdii,fDMi], MDM];
AppendTo[PtDMList[c\[Nu]dii,csdii,fDMi],PteDM[\[Rho]DM,MDM,kDM]];
AppendTo[\[Epsilon]tDMList[c\[Nu]dii,csdii,fDMi],\[Epsilon]teDM[\[Rho]DM,MDM,kDM]];
]
];
];
];(*EndFor*)
(*############################################################################################################################*)
rhodata=\[Rho]List[c\[Nu]d,csd,fDMi](*[[;;-50]]*);
\[Epsilon]NMdata=\[Epsilon]tList[c\[Nu]d,csd,fDMi](*[[;;-50]]*);
\[Epsilon]DMdata=\[Epsilon]tDMList[c\[Nu]d,csd,fDMi](*[[;;-50]]*);
PNMdata=PtList[c\[Nu]d,csd,fDMi](*[[;;-50]]*);
PDMdata=PtDMList[c\[Nu]d,csd,fDMi](*[[;;-50]]*);

Final=Transpose[
{\[Rho]List[c\[Nu]dii,csdii,fDMi],
\[Epsilon]tList[c\[Nu]dii,csdii,fDMi],
PtList[c\[Nu]dii,csdii,fDMi],
\[Epsilon]tDMList[c\[Nu]dii,csdii,fDMi],
PtDMList[c\[Nu]dii,csdii,fDMi]}];
(*SetDirectory["/home/everson/Documents/DM with SRC/Data-SRC/Bosonic-Md=15000-Data"];

file1=OpenWrite["PtandEt-csd="<>ToString[NumberForm[csdii,{3,2}]]<>", fr="<>ToString[NumberForm[fDMi,{5,4}]]<>", SR="<>ToString[SRi]<>", "<>modelname<>"-test.dat",BinaryFormat->True];
Export[file1,Final];
Close[file1];*)

(*
CompoundExpression[]
Byte count: 9176
Uniconize


*)
Return[Final]];
End[]
EndPackage[]

