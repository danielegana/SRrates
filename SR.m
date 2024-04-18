(* ::Package:: *)

(* ::Input::Initialization:: *)
directory="/home/degana/SR/";
SetDirectory[directory];


(* ::Section::Closed:: *)
(*Preliminaries*)


(* ::Subsection::Closed:: *)
(*Conversions*)


(* ::Input::Initialization:: *)
barntocm2=10^(-28)*10^4;
cm2tobarn=1/barntocm2;
cm2topbarn=cm2tobarn*10^12;
eVtokg=1.78*10^(-36);
GeVtokg=1.78*10^(-27);

kgtoGeV=1/GeVtokg;
gtoGeV=kgtoGeV/1000;

invGeVtocm=1.97`*^-14;
invGeVtom=1.97`*^-16;
inveVtocm=1.97`*^-5;
inveVtokm=1.97`*^-10;
invGeVtoMPc=invGeVtocm/(100*3.086*10^22);
invGeVtokPc=invGeVtoMPc*10^3;
invGeVtoPc=invGeVtoMPc*10^6;
GeVtoinvMPc=1/invGeVtoMPc;
MPctoinvGeV=1/invGeVtoMPc;
PctoinvGeV=MPctoinvGeV*10^(-6);
Pctom=30.9*10^15;
Pctocm=30.9*10^17;
kPctoinvGeV=MPctoinvGeV*10^(-3);
invMPctoGeV=1/MPctoinvGeV;

kmtoinveV=1/inveVtokm;
kmtoinvGeV=10^9*kmtoinveV;
cmtoinvGeV=1/invGeVtocm;
invGeVtokm=1/kmtoinvGeV;
nmtoinvGeV=10^(-7)*cmtoinvGeV;
\[Mu]mtoinvGeV=10^(-4)*cmtoinvGeV;

invcmtoGeV=1/cmtoinvGeV;
GeVtoinvcm=1/invcmtoGeV;

eVtokelvin=11600;
GeVtokelvin=10^9*11600;
kelvintoeV=1/11600;
kelvintoGeV=kelvintoeV/10^9;

inveVtos=6.58*10^(-16);
invGeVtos=6.58*10^(-16)*10^(-9);
GeVtoinvs=1/invGeVtos;
eVtoHz=10^(-9)*GeVtoinvs;
eVtoMHz=10^(-15)*GeVtoinvs;

stoinveV=1/inveVtos;
stoinvGeV=1/invGeVtos;
ytoinveV=3.154*10^7*stoinveV;
ytoinvGeV=ytoinveV*10^9;
GytoinveV=10^9*ytoinveV;
GytoinvGeV=10^9*GytoinveV;

invstoeV=1/stoinveV;
invstoGeV=1/stoinvGeV;
invGytoGeV=1/GytoinvGeV;
invGeVtoGy=1/GytoinvGeV;
invGeVtoy=invGeVtoGy*10^9;

ergtoGeV=624.15;
GeVtoerg=1/ergtoGeV;
stoyear=1/(60*60*24*365);
stoday=1/(60*60*24);
Coulomb=6.2*10^(18);
TeslatoGeV2=kgtoGeV*invstoGeV/(Coulomb);


(* ::Subsection:: *)
(*Numbers and definitions*)


(* ::Subsubsection:: *)
(*Particle physics*)


(* ::Input::Initialization:: *)
replnumbers=Rationalize[{\[Alpha]fine->1/137.035999074,e->Sqrt[1/137.035999074*4*Pi],\[Theta]w->ArcSin[Sqrt[0.23126]],mW->80.376,mZ->91.185,cw->0.8767781931594786`,sw->Sqrt[0.23126],v->246,g->0.6491458936073573`,gp->0.35604331431828873`,s->207.5^2,\[CapitalGamma]Z->2.5,ed->-1,me->0.511*(10)^(-3),m\[Pi]->0.1396,Gf->1.16*10^(-5),f\[Pi]->0.130,\[Theta]c->ArcSin[0.22],MPl->2.43*10^(18),mPl->Sqrt[8*Pi]*2.43*10^(18),G->6.738257341017432`*^-39(*In GeV^-2*),Mmilkiway->10^12*MSolGeV,H0invGY->0.673/9.77},0];
RnucleoninvGeV=1.25*cmtoinvGeV*10^(-13);
\[Rho]nucl=1/(4*Pi*RnucleoninvGeV^3/3);
\[Rho]nuclkgcm3=1/(4*Pi*RnucleoninvGeV^3/3)*GeVtokg/invGeVtocm^3;
EFnucl=1/2*(3*Pi^2*\[Rho]nucl)^(3/2);
EFNS=1/2*(3*Pi^2*2*\[Rho]nucl)^(3/2);
pFnucl=Sqrt[2*EFnucl];
pFNS=Sqrt[2*EFNS];


(* ::Subsubsection:: *)
(*Astro parameters*)


(* ::Input::Initialization:: *)
MSolkg=1.99*10^(30);
RSuninvGeV=695508*kmtoinvGeV;
MSolGeV=Rationalize[MSolkg*kgtoGeV,0];
MmilkiwayGeV=MSolGeV*10^12;
Csun=G*MSolGeV/RSuninvGeV/.replnumbers;
distancesuncm=149.6*10^6*10^3*10^2;

RSS[MBHSol_]:=2*G*MBHSol*MSolGeV/.replnumbers;
sqdgsky=360^2/Pi;


(* ::Subsubsection:: *)
(*Cosmological parameters*)


(* ::Input::Initialization:: *)

(*FIRAS, Fixsen 0911.1955*)
TCMBkelvin=Rationalize[2.7255,0];
TCMBGeV=Rationalize[TCMBkelvin*kelvintoeV/10^9,0];
s0=2*Pi^2/45*3.91*TCMBGeV^3;

g\[Gamma]=2;g\[Nu]=6;T\[Nu]=TCMBGeV*(4/11)^(1/3);
e\[Gamma][T\[Gamma]_]:=g\[Gamma]*Pi^2/30*(T\[Gamma])^4;
e\[Nu][T_]:=g\[Nu]*7/8*Pi^2/30*(T)^4; (*=Neff*7/8*(4/11)^(4/3)*e\[Gamma]*)


(*Planck 2018+BAO+Lensing*)
hcosmo=0.6766; 
H0GeV=100*kmtoinvGeV/stoinvGeV*invMPctoGeV*hcosmo;
\[Rho]c=3*H0GeV^2/(8*Pi*G)/.replnumbers;

\[CapitalOmega]\[Gamma]o=e\[Gamma][TCMBGeV]/\[Rho]c;

\[CapitalOmega]ro=4.2*10^(-5)/hcosmo^2; (*Not updated with new CMB temp*)
\[CapitalOmega]mo=0.11933/hcosmo^2+0.02242/hcosmo^2;
\[CapitalOmega]\[CapitalLambda]o=1-\[CapitalOmega]mo-\[CapitalOmega]ro;

\[Rho]DM=Rationalize[0.11933/hcosmo^2*\[Rho]c,0];
\[Rho]b=Rationalize[0.02242/hcosmo^2*\[Rho]c,0];



(*Approximate numbers *)
zeq=2740;zrec=1100;
\[Rho]DMlocalGeVovercm3=Rationalize[0.3,0];
\[Rho]DMlocal=Rationalize[0.3/cmtoinvGeV^3,0];

(*Hubble*)
Efun[z_]:=Sqrt[\[CapitalOmega]\[CapitalLambda]o+(1-\[CapitalOmega]\[CapitalLambda]o-\[CapitalOmega]mo-\[CapitalOmega]ro)*(1+z)^2+\[CapitalOmega]mo*(1+z)^3+\[CapitalOmega]ro*(1+z)^4]
HzGeV[z_]:=H0GeV*Efun[z]; (*H^2=8PiG/3*\[Rho](z). H0^2=8*Pi*G/3 \[Rho]c. So one can either choose rhoc, H0 or h as an independent cosmo parameter.*)
HzradGeV[z_]:=H0symb*Sqrt[\[CapitalOmega]rosymb]*(1+z)^2;
HaGeV[a_]:=HzGeV[(1-a)/a];
H\[Eta]GeV[a_]:=HzGeV[(1-a)/a]*a;
HaradGeV[a_]:=(H0symb Sqrt[\[CapitalOmega]rosymb])/a^2;

(*Densities*)
\[Rho]GeV[z_]:=HzGeV[z]^2*3/(8*Pi*G);
\[Rho]aGeV[a_]:=HaGeV[a]^2*3/(8*Pi*G);
\[Rho]DMa[a_]:=\[Rho]DM/a^3;
\[Rho]DMz[z_]:=\[Rho]DMa[1/(1+z)];
\[Rho]ba[a_]:=\[Rho]b/a^3;
\[Rho]Ma[a_]:=(\[Rho]DM+\[Rho]b)/a^3;
\[Rho]rada[a_]:=\[CapitalOmega]ro*\[Rho]c/a^4;
azfun[z_]:=1/(1+z);
zafun[a_]:=(1-a)/a;

\[CapitalOmega]wdm[mkeV_,TdecoverT_,g_]:=15*mkeV*TdecoverT^3*g/hcosmo^2;


(* ::Text:: *)
(*Matter DE equality*)


(* ::Input::Initialization:: *)
replnumcosmo={H0symb->H0GeV,\[CapitalOmega]rosymb->\[CapitalOmega]ro};


(* ::Subsubsection:: *)
(*Thermodynamics *)


(* ::Input::Initialization:: *)
\[Sigma]SB=2*Pi^5/(15*(2*Pi)^3);
\[Nu]max[T_]=2.821439372122078893*T/(2*Pi);
piB\[Nu][\[Nu]_,T_]:=Rationalize[Pi*2*(2*Pi),0]*\[Nu]^3/(Exp[2*Rationalize[Pi,0]*\[Nu]/T]-1);(*Factor of Pi in front from solid angle integral*)
IBB[T_?NumericQ]:=Rationalize[\[Sigma]SB,0]*T^4;


(* ::Section:: *)
(*Packages and data import*)


(* ::Subsection:: *)
(*Find all roots*)


(* ::Input::Initialization:: *)
Clear[findAllRoots]
SyntaxInformation[findAllRoots]={"LocalVariables"->{"Plot",{2,2}},"ArgumentsPattern"->{_,_,OptionsPattern[]}};
SetAttributes[findAllRoots,HoldAll];

Options[findAllRoots]=Join[{"ShowPlot"->False,PlotRange->All},FilterRules[Options[Plot],Except[PlotRange]]];

findAllRoots[fn_,{l_,lmin_,lmax_},opts:OptionsPattern[]]:=Module[{pl,p,x,localFunction,brackets},localFunction=ReleaseHold[Hold[fn]/.HoldPattern[l]:>x];
If[lmin!=lmax,pl=Plot[localFunction,{x,lmin,lmax},Evaluate@FilterRules[Join[{opts},Options[findAllRoots]],Options[Plot]]];
p=Cases[pl,Line[{x__}]:>x,Infinity];
If[OptionValue["ShowPlot"],Print[Show[pl,PlotLabel->"Finding roots for this function",ImageSize->200,BaseStyle->{FontSize->8}]]],p={}];
brackets=Map[First,Select[(*This Split trick pretends that two points on the curve are "equal" if the function values have _opposite _ sign.Pairs of such sign-changes form the brackets for the subsequent FindRoot*)Split[p,Sign[Last[#2]]==-Sign[Last[#1]]&],Length[#1]==2&],{2}];
x/.Apply[FindRoot[localFunction==0,{x,##1}]&,brackets,{1}]/.x->{}]


(* ::Section:: *)
(*Accretion rates (Run!)*)


(* ::Input::Initialization:: *)
\[Sigma]T=8*Pi/3*(\[Alpha]fine/me)^2/.replnumbers;
Ledd[Msol_]=4*Pi*G*Msol*MSolGeV/\[Sigma]T;
Lacc[\[Eta]_]:=\[Eta]*Mdot;
Mdotedd[Msol_,\[Eta]_]:=(8.223098154932542`*^54 G Msol)/\[Eta];
Mdoteddsolyear[Msol_,\[Eta]_]:=Mdotedd[Msol,\[Eta]]*GeVtoinvs/stoyear/MSolGeV/.replnumbers;
\[Tau]eddyear=(Msol*MSolGeV)/Ledd[Msol]*invGeVtoGy*10^9/.replnumbers
4.2093419872719586`*^8


(* ::Section::Closed:: *)
(*Spin up equations (run)*)


(* ::Subsection:: *)
(*Isco *)


(* ::Input::Initialization:: *)
Z1[a_]:=1+(1-a^2)^(1/3)*((1+a)^(1/3)+(1-a)^(1/3));
Z2[a_]:=(3*a^2+Z1[a]^2)^(1/2);
riscooverrg[a_]:=(3+Z2[a]-Sqrt[(3-Z1[a])*(3+Z1[a]+2*Z2[a])]);

liscoovergM[a_]:=Sqrt[riscooverrg[a]]*(riscooverrg[a]^2-2*a*Sqrt[riscooverrg[a]]+a^2)/(riscooverrg[a]*Sqrt[riscooverrg[a]^2-3*riscooverrg[a]+2*a*Sqrt[riscooverrg[a]]]);


(* ::Subsection:: *)
(*Accretion efficiency*)


(* ::Input::Initialization:: *)
Ensol[adimless_]:=Module[{Entemp},
Entemp=(En/.Solve[adimless==(4*Sqrt[2]*(1-En^2)^(1/2)-2*En)/(3*Sqrt[3]*(1-En^2)),En]);
Select[Entemp,#\[Element]Reals&][[1]]
]
\[Epsilon][adimless_]:=1-Ensol[adimless];
maxspinfactor[adimless_]:=2*(Sqrt[1-adimless^2]+1-adimless^2)/adimless;
sover\[Epsilon]areatheorem[adimless_]:=2*(1-\[Epsilon][adimless])/\[Epsilon][adimless]*(Sqrt[1-adimless^2]+1-adimless^2)/adimless;
sover\[Epsilon][adimless_]:=If[adimless==1,0,(1-\[Epsilon][adimless])/\[Epsilon][adimless]*(liscoovergM[adimless]/(1-\[Epsilon][adimless])-2*adimless)];
atable=Range[0,1,0.001];
sover\[Epsilon]int=Interpolation[Transpose[{atable,ParallelMap[sover\[Epsilon][#]&,atable]}]];


(* ::Section:: *)
(*Superradiance rates*)


(* ::Input::Initialization:: *)
rg[MBHSol_]:=G*MBHSol*MSolGeV;
rss[MBHSol_]:=2*G*MBHSol*MSolGeV;
\[CapitalOmega]kerrada[MBHSol_,a_]:=1/2*(a/(G*MBHSol*MSolGeV)/(1+Sqrt[1-(a/(G*MBHSol*MSolGeV))^2]))/rg[MBHSol];(*Equal to a/(2*G*M*rkerr)*)
\[CapitalOmega]kerr[MBHSol_,adimless_]:=1/2*(adimless/(1+Sqrt[1-adimless^2]))/rg[MBHSol];
\[Alpha][\[Mu]_,MBHSol_]:=\[Mu]*rg[MBHSol]/.replnumbers;
\[Mu]\[Alpha][\[Alpha]_,MBHSol_]:=\[Alpha]/rg[MBHSol]/.replnumbers;
M\[Alpha]sol[\[Alpha]_,\[Mu]_]:=\[Alpha]/(\[Mu]*G)/MSolGeV/.replnumbers;
rkerr[MBHsol_,a_]:=(rss[MBHsol]+Sqrt[rss[MBHsol]^2-4*a^2])/2;

\[CapitalGamma][\[Mu]_,MBHsol_,n_,l_?NumericQ,m_,a_]:=\[Mu]*(\[Mu]*G*MBHsol*MSolGeV)^(4*l+4)*(a*m/(G*MBHsol*MSolGeV)-2*\[Mu]*rkerr[MBHsol,a])*2^(4*l+2)*((2*l+1+n)!)/(((l+1+n)^(2*l+4))*n!)*
(l!/(((2*l)!)*((2*l+1)!)))^2*
\!\(
\*SubsuperscriptBox[\(\[Product]\), \(j = 1\), \(l\)]\((j^2*\((1 - a^2/\((G*MBHsol*MSolGeV)\)^2)\) + \[IndentingNewLine]\((a*m/\((G*MBHsol*MSolGeV)\) - 2*\[Mu]*rkerr[MBHsol, a])\)^2)\)\)/.replnumbers; (*n starts from 0. n is not the principal quantum number, it's what Shankar calls k in page 355. The im part of the frequency \[Omega] is \[CapitalGamma]/2.*)

\[CapitalGamma]sm[\[Mu]_,MBHsol_,n_,l_,m_,a_]:=1/2*
\[Mu]*(\[Mu]*G*MBHsol*MSolGeV)^(4*l+4)*(2*m*\[CapitalOmega]kerrada[MBHsol,a]*rkerr[MBHsol,a]-2*\[Mu]*rkerr[MBHsol,a])*2^(4*l+2)*((2*l+1+n)!)/(((l+1+n)^(2*l+4))*n!)*
(l!/(((2*l)!)*((2*l+1)!)))^2*
\!\(
\*SubsuperscriptBox[\(\[Product]\), \(j = 1\), \(l\)]\((j^2*\((1 - a^2/\((G*MBHsol*MSolGeV)\)^2)\) + \[IndentingNewLine]\((2*m*\[CapitalOmega]kerrada[MBHsol, a]*rkerr[MBHsol, a] - 2*\[Mu]*rkerr[MBHsol, a])\)^2)\)\)/.replnumbers(*Factor of 1/2 missing in original reference*)

\[CapitalGamma]adimless[\[Mu]_,MBHsol_,n_,l_,m_,adimless_]:=\[CapitalGamma]sm[\[Mu],MBHsol,n,l,m,adimless*(G*MBHsol*MSolGeV)];
SRcondition\[Mu][MBHsol_,m_,adimless_]:=m*\[CapitalOmega]kerr[MBHsol,adimless]/.replnumbers;
logNmax[MBHsol_,m_,da_]:=Log10[G*(MBHsol*MSolGeV)^2/m*da/.replnumbers];

asolSRtemp=adimless/.Solve[\[Mu]p==mp*\[CapitalOmega]kerr[MBHsolp,adimless],adimless][[1]];
SRconditiona[MBHsol_,m_,\[Mu]_]:=asolSRtemp/.{MBHsolp->MBHsol,mp->m,\[Mu]p->\[Mu]}/.replnumbers;
\[Tau]SRyears[\[Mu]_,MBHsol_,n_,l_,m_,adimless_]:=If[\[Mu]<SRcondition\[Mu][MBHsol,m,adimless],(\[CapitalGamma]adimless[\[Mu],MBHsol,n,l,m,adimless]*GeVtoinvs/stoyear)^(-1),\[Infinity]];

\[Tau]ref[adimless_?NumericQ]:=If[adimless==1,\[Infinity],\[Tau]eddyear/sover\[Epsilon]int[adimless]];

magneticmax=4;

ageuniverseyears=1.3813246371323896`*^10;
(*Compare with disk spinup time. Only take n=0.*)
(*aextr[\[Mu]_?NumericQ,MBHsol_?NumericQ,m_?NumericQ]:=
Module[{\[Mu]0=\[Mu],MBHsol0=MBHsol,m0=m,aloc},
If[\[Mu]>SRcondition\[Mu][MBHsol,m,1]||\[Tau]SRyears[\[Mu],MBHsol,0,m,m,1]>ageuniverseyears,aloc=1,
aloc=adimless/.FindRoot[logNmax[MBHsol,m,0.1]*\[Tau]SRyears[\[Mu],MBHsol,0,m,m,adimless]==\[Tau]ref[adimless],{adimless,0.999,1.00001*SRconditiona[MBHsol,m,\[Mu]],1}][[1]]];
(*aloc=adimless/.FindRoot[logNmax[MBHsol,m,0.1]*\[Tau]SRyears[\[Mu],MBHsol,m+1,m,m,adimless]-\[Tau]ref,{adimless,0}][[1]];*)
If[aloc>0,aloc,1,1]
];*)

(*Compare with disk spinup time. Take n=0 and n=1.*)
aextr[\[Mu]_?NumericQ,MBHsol_?NumericQ,m_?NumericQ]:=
Module[{\[Mu]0=\[Mu],MBHsol0=MBHsol,m0=m,aloc},
aloc=Min[adimless/.FindRoot[logNmax[MBHsol,m,0.1]*\[Tau]SRyears[\[Mu],MBHsol,0,m,m,adimless]==\[Tau]ref[adimless],{adimless,0.999,1.001*SRconditiona[MBHsol,m,\[Mu]],1}][[1]],adimless/.FindRoot[logNmax[MBHsol,m,0.1]*\[Tau]SRyears[\[Mu],MBHsol,1,m,m,adimless]==\[Tau]ref[adimless],{adimless,0.999,1.001*SRconditiona[MBHsol,m,\[Mu]],1}][[1]]];
(*aloc=adimless/.FindRoot[logNmax[MBHsol,m,0.1]*\[Tau]SRyears[\[Mu],MBHsol,m+1,m,m,adimless]-\[Tau]ref,{adimless,0}][[1]];*)
If[aloc>0,aloc,1,1]
];
aexttable[\[Mu]_,MBHsol_]:=Quiet[Table[{m,Quiet[aextr[\[Mu],MBHsol,m]//N]},{m,1,magneticmax}]];
aextmin[\[Mu]_?NumericQ,MBHsol_?NumericQ]:=SortBy[aexttable[\[Mu],MBHsol],#[[2]]&][[1]];




(* ::Section:: *)
(*Hills table*)


(* ::Input::Initialization:: *)
(*atable=Flatten[{Range[0,0.9,0.1],{0.95,0.99}}];
logMHilltable=ParallelMap[logMHillssol\[Iota][#,0,1,RSuninvGeV]&,atable];
MHillint=Interpolation[Transpose[{atable,10^logMHilltable}],InterpolationOrder\[Rule]1];
aHillint=Interpolation[Transpose[{10^logMHilltable,atable}],InterpolationOrder\[Rule]1];*)
logMHilltable=<<"logMHilltable.txt";


(* ::Input:: *)
(*(*logMHilltable<<"logMHilltable.txt"*)*)


(* ::Section:: *)
(*Axion SR tables*)


(* ::Input::Initialization:: *)
\[Mu]min=10^(-21)*10^(-9);\[Mu]max=10^(-18)*10^(-9);
nclus=
ntot=20;
nmult=2000;
nben=ntot*nmult;
nbentable=Range[1,nben];
log\[Mu]bentable=Range[Log10[\[Mu]min],Log10[\[Mu]max],(Log10[\[Mu]max]-Log10[\[Mu]min])/(nben-1)][[(nclus-1)*nmult+1;;nclus*nmult]];
\[Mu]bentable=Rationalize[(10^#)&/@log\[Mu]bentable,0];
dim\[Mu]=Dimensions[\[Mu]bentable][[1]];

tablespin=Table[0,{n,1,dim\[Mu]}];

mmin=10^6;mmax=10^9;
npoints=2000;
logntable=N[Range[Log10[mmin],Log10[mmax],(Log10[mmax]-Log10[mmin])/npoints]];
ntable=Rationalize[(10^#)&/@logntable,0];

Do[tablespin[[n]]=ParallelMap[aextmin[\[Mu]bentable[[n]],#]&,ntable];If[n==1,Print[tablespin[[n]]],Print[n]];,{n,1,dim\[Mu]}];
Put[tablespin,StringJoin["output/tablespin",ToString[nclus],".dat"]];
Put[\[Mu]bentable,StringJoin["output/\[Mu]bentable",ToString[nclus],".dat"]];
