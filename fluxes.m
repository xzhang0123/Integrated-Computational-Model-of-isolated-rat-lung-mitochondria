function rtfluxes=fluxes(C_t,Para1)


global k0 Tem F_con R_con Ve Vm Vi ROTi AAi closed_system...
 iPie  iADPe iATPe iPYRe iMALe iCITe iaKGe iSUCe iFUMe iGLUe iASPe...
    iGLUm iASPm iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iNADm iNADHm...
    iUQm iUQH2m iCytCoxi iCytCred iADPm iATPm iGDPm iGTPm iCOAm iACOAm iPim...
    iFADm iFADH2m iHm iHe idPsi iO2 iR123e iR123m iReBe iRmBm
RT=R_con*Tem;

dGr(1)= -38.64; %PYR dehydrogenase [1]
dGr(2)= -36.61;%CIT synthase [1]
dGr(3)=2.18;%aconiase+isocit dehydrogenase [1]
dGr(4)=-37.08;%aKG dehydro
dGr(5)=1.26;%scoa synthetase
dGr(6)=0;   % GDP+ATP-GTP-ADP
dGr(7)=-3.62;%suc dehydro
dGr(8)=-3.6; %FUM Jason's model (not generating H+, so assumed not changing with pH)
dGr(9)=28.83;%mal dehydro
dGr(10)=-1.4;%asp+akg-glu+oxa [3]
dGr(11)=-69.37;%CI   [2] 
dGr(12)= -1.31;   %[2]-[1]
dGr(13)=-32.53; %CIII  [2]
dGr(14)=-122.94;%CIV   [2]
dGr(15)=36.03;%ATPase [2]
%referneces:
%[1] Role of NADH/NAD transport 
%[2] analysis of cardiac mito NA/Ca exchange kinetics with a biophysical
%model of ~~~~ 3:1 stoichiometry
%[3] Appendices for ?Computer Modeling of Mitochondrial TCA Cycle, Oxidative
%Phosphorylation, Metabolite Transport, and Electrophysiology?
Keq0=exp(-dGr/RT);

  Est_Para=ones(25,1);
 mito_yield_ratio=1300/2;
Est_Para=[306.49e-09;
    4630.8*1e-9;
    530.797548526789e-09;
    %3.29889721588679e3*1e-9;% AKGDH
    3*mito_yield_ratio*0.4e-9;% AKGDH ref
    5.81695925037748e3*1e-9;
    4.32513873910319e5*1e-09;  %NDK not identifiable
    1.08593076185781e-05;   %SDH
    mito_yield_ratio*9.9e-9; %FH ref=Sources of C for TCA 
    3297.559540729090*1e-9;    %MDH
   	mito_yield_ratio* 1.5e-9; %got,ref=Sources of C for TCA 
    11.007e-9;      %CI
   220*1e-9;      %CII  ref=Age related changes in metabolism oxidative stress
 2.2335e4*1e-9;        %CIII
 0.271043306907577*1e-9;        %CIV
    588.98e-09;    %CV
    %----------------------------------
    1698.5e-09;  %Tmax1 DCC SUC
    21.98e-09;    %DCC MAL
    0.9994e-09;   %OME
   1/4* 81.33e-09;   %TCC
    96.62e-09;      %PYRH
    2.37958063427568e4*1e-9;     %PIC
   523.8837274961545e-09;      %ANT
    4.626608507608435e-9;        %GLUH
   646.25*1e-9;       %asp-glu
  35.5*1e-9;
    0;      %Not used
    0];

% Vmaxf=Para1(1:15);
% Tmax=Para1(16:27);
Vmaxf=Para1(1:15).*Est_Para(1:15);
Tmax=Para1(16:25).*Est_Para(16:25);

CO2=1.32e-5;   %henry's law
%-----
% % M(28:29)=4e-3;
% % M(30:31)=0.7e-3;
% % M(32:33)=1e-4;
% % M(34)=2e-3;
%-----------
Pie  = C_t(1);
ADPe = C_t(2);
ATPe = C_t(3);
PYRe = C_t(4);
MALe = C_t(5);
CITe = C_t(6);
aKGe = C_t(7);
SUCe = C_t(8);
FUMe = C_t(9);
GLUe = C_t(10);
ASPe = C_t(11);
GLUm = C_t(12);
ASPm = C_t(13);
PYRm = C_t(14);
OXAm = C_t(15);
CITm = C_t(16);
aKGm = C_t(17);
SCAm = C_t(18);
SUCm = C_t(19);
FUMm = C_t(20);
MALm = C_t(21);
NADm  = C_t(22);
NADHm = C_t(23);
UQm   = C_t(24);
UQH2m =C_t(25);
CytCoxi = C_t(26);
CytCred = C_t(27);
ADPm   = C_t(28);
ATPm = C_t(29);
GDPm = C_t(30);
GTPm = C_t(31);
COAm= C_t(32);
ACOAm= C_t(33);
Pim = C_t(34);
FADm =C_t(35);
FADH2m = C_t(36);    
Hm=C_t(37);
He=C_t(38);
dPsi=C_t(39);
O2=C_t(40);
R123e=C_t(41);
R123m=C_t(42);
ReBe=C_t(43);
RmBm=C_t(44);
DeltaGH = (F_con*dPsi+ R_con*Tem*log(He/Hm));

pH_e=-log10(He);
pH_m=-log10(Hm);
%----------------------------------------------
%R1: PYR+COAm+NAD->ACOA+CO2+NADH+Hm
%KA1=2.52e-5; %change to estiamted paremeter
KA1=5.08e-3;
KB1=1.49e-5;
KC1=3.5e-5;
KD1=1.49e-5;
%KD1=8e-5;
KE1=1e-5;%not used in model
KF1=3.5e-5;
A=PYRm;
B=COAm;
C=NADm;
D=ACOAm;
E=CO2;
F=NADHm;
Keq(1)=Keq0(1)*10^(pH_m-7);
deno=(1+A/KA1)*(1+B/KB1+D/KD1)*(1+C/KC1+F/KF1);   %(1+1+E/KE) is ignored since it is constant
%--------------------------
%---------------------------
J_PDH =Vmaxf(1)/KA1/KB1/KC1*(A*B*C-D*E*F/Keq(1))/deno;
%-----------------------------------------------
%R2: ACOA+OXA->COA+CIT+Hm
% KA2=4.6e-6;
% KB2=5.7e-6;
% KC2=3.5e-3;
% KD2=4.9e-5;
KA2=3.9e-6;
KB2=4.53e-6;
KC2=57.9e-6;
KD2=4.3e-3;
A=ACOAm;
B=OXAm;
C=COAm;
D=CITm;
Keq(2)=Keq0(2)*10^(pH_m-7)*10^(pH_m-7);
deno=(1+A/KA2+C/KC2)*(1+B/KB2)*(1+D/KD2);
J_CITS =Vmaxf(2)/KA2/KB2*(A*B-C*D/Keq(2))/deno;
%----------------------------------------------
%R3: CIT+NAD+->aKGm+NADHm+CO3--+Hm
KA3=1.3e-3;
KB3=500e-6;
KC3=3.5e-6;
KD3=4.7e-6;
KE3=1e-5; %CO2,NOT USED
%--------------------
% KA3=1.3e-6;
% KB3=500e-6;
% KC3=3.5e-6;
% KD3=49e-6;
% KE3=1e-5; %CO2,NOT USED
%-----------
A=CITm;
B=NADm;
C=aKGm;
D=NADHm;
E=CO2;
Keq(3)=Keq0(3);
deno=(1+A/KA3)*(1+C/KC3)*(1+B/KB3+D/KD3);
J_CITD =Vmaxf(3)/KA3/KB3*(A*B-C*D*E/Keq(3))/deno;
%----------------------------------------------
%R4: aKG + CoA+NAD+H2O-> SCoA+ NADH+H2CO3+H+
% KA4=120e-6;    %AKG
% KB4=4e-7;      %COA
% KC4=29e-6;     %NAD
% % KD4=3.0e-3;    %SCA
% KD4=KB4;    %SCA
% KE4=0.15e-6;   %NADH 
% %KE4=KC4;   %NADH 
% KF4=1e-5;%CO2,NOT USED, co2 constant
%-------
KA4=85.87e-6;    %AKG
KB4=1.634e-6;      %COA
KC4=46.6e-6;     %NAD
% KD4=3.0e-3;    %SCA
KD4=KB4;    %SCA
KE4=7.14e-6;   %NADH 
%KE4=KC4;   %NADH 
KF4=1e-5;%CO2,NOT USED
%--------------------
A=aKGm;
B=COAm;
C=NADm;
D=SCAm;
E=NADHm;
F=CO2;
Keq(4)=Keq0(4)*10^(7-pH_m);
deno=(1+B/KB4+D/KD4)*(1+C/KC4+E/KE4)*(1+A/KA4);
J_AKGD =Vmaxf(4)/KA4/KB4/KC4*(A*B*C-D*E*F/Keq(4))/deno;
%----------------------------------------------
%R5 SCA+GDPm+Pim ->SUC+GTP+COA+H+
KA5=1.2538e-5;
KB5=5.6977e-6;
KC5=2.1e-3;
KD5=4.51e-4;
KE5=2.3143e-5;
KF5=1.7735e-5;
A=SCAm;
B=GDPm;
C=Pim;
D=SUCm;
E=GTPm;
F=COAm;
Keq(5)=Keq0(5)*10^(pH_m-7);
%deno=(1+A/KA5+D/KD5+E/KE5+D*E/KD5/KE5)*(1+B/KB5+C/KC5+F/KF5+B*C/KB5/KC5);
deno=(1+A/KA5+D/KD5+E/KE5+D*E/KD5/KE5)*(1+B/KB5+C/KC5+F/KF5+B*C/KB5/KC5);
J_SCAS= Vmaxf(5)/KA5/KB5/KC5*(A*B*C-D*E*F/Keq(5))/deno;
%----------------------------------------------
%R6 GTP+ADP ----GDP+ATP
KA6=111e-6;
KB6=100e-6;
KC6=260e-6;
KD6=278e-6;
A=GTPm;
B=ADPm;
C=GDPm;
D=ATPm;
Keq(6)=Keq0(6);
deno=(1+A/KA6+C/KC6)*(1+B/KB6+D/KD6);
J_NDK=Vmaxf(6)/KA6/KB6*(A*B-C*D/Keq(6))/deno;
%----------------------------------------------
%R7 SUC+FAD->FUM+FADH2
KA7=1800e-6;
KB7=140e-6;
KC7=1800e-6;
KD7=2.45e-6;
A=SUCm;
B=FADm;
C=MALm;
D=FADH2m;
Keq(7)=Keq0(7);
deno=(1+A/KA7)*(1+C/KC7)*(1+B/KB7+D/KD7);
J_SUCD =Vmaxf(7)/KA7/KB7*(A*B-C*D/Keq(7))/deno;
%-----------------------------------------------
%R8 FUM+H2O->MAL
KA8=1800e-6;
KB8=1800e-6;

A=FUMm;
B=MALm;

Keq(8)=Keq0(8);
deno=(1+A/KA8)*(1+B/KB8);
J_FUM =Vmaxf(8)/KA8/KB8*(A-B/Keq(8))/deno;
%----------------------------------------------
%R9 MAL+NAD+->OXA+NADH+Hm

KA9=1.55e-4;
KB9=1.1e-3;
KC9=3.38e-6;
KD9=3.61e-5;
A=MALm;
B=NADm;
C=OXAm;
D=NADHm;
Keq(9)=Keq0(9)*10^(pH_m-7);
deno=(1+A/KA9+C/KC9)*(1+B/KB9+D/KD9);
% deno=(1+A/KA9)*(1+C/KC9)*(1+B/KB9+D/KD9);
J_MALD =Vmaxf(9)/KA9/KB9*(A*B-C*D/Keq(9))/deno;
%----------------------------------------------
% R10 ASPm+aKGm-GLUm+OXA
KA10=3.9e-3;  %use jason's parameters
KB10=430e-6;   %AKG
KC10=8.9e-3;   %GLU
KD10=88e-6;     %OAA

A=ASPm;
B=aKGm;
C=GLUm;
D=OXAm;
Keq(10)=Keq0(10);
deno=(1+A/KA10)*(1+B/KB10)*(1+C/KC10)*(1+D/KD10);
J_GOT =Vmaxf(10)/KA10/KB10*(A*B-C*D/Keq(10))/deno;
% %----------------------------------------------
 betaC1=0.5;
%-----------------------
% betaC1=0;
% Vmaxf(11)=Vmaxf(11)*200e-6;
%----------------------
KA11=1.5e-6;
KB11=58.1e-6;
KC11=428e-6;
KD11=520e-6;
A=NADHm;
B=UQm;
C=NADm;
D=UQH2m;
Keq(11)=Keq0(11)*10^(7-pH_m);
deno=(1+A/KA11+C/KC11)*(1+B/KB11+D/KD11);
%5e4-1e5 works
% exp(-2*DeltaGH/RT)
a=Keq(11).^0.5;
J_CI =ROTi*Vmaxf(11)/KA11/KB11*(a*exp(-4*betaC1*DeltaGH/RT)*A*B-exp(-4*(betaC1-1)*DeltaGH/RT)*C*D/a)/deno;
% %----------------------------------------------
% %R12 Complex II    FADH2+UQm--FADm+UQH2m
KA12=1.5e-6;
KB12=58e-6;
KC12=428e-6;
KD12=519e-6;
% KA12=24.2e-6;
% KB12=58e-6;
% KC12=0.12e-6;
% KD12=520e-6;
A=FADH2m;
B=UQm;
C=FADm;
D=UQH2m;
Keq(12)=Keq0(12);
deno=(1+A/KA12+C/KC12)*(1+B/KB12+D/KD12);
J_CII=Vmaxf(12)/KA12/KB12*(A*B-C*D/Keq(12))/deno;
% %----------------------------------------------
% %R13 Complex III   UQH2m+2CytCoxi+2Hm--UQm+2CytCred+4Hi
 betaC3=0.5;
%------------
KA13=4.66e-6;
KB13=3.76e-6;
KC13=4.08e-6;
KD13=4.91e-6;
A=UQH2m;
B=CytCoxi;
C=UQm;
D=CytCred;
Keq(13)=Keq0(13)*10^(pH_m-7)*10^(pH_m-7);%
%whole reaction
deno=(1+A/KA13+C/KC13)*(1+B^2/KB13^2+D^2/KD13^2);
J_CIII=AAi*Vmaxf(13)/KA13/(KB13^2)*(Keq(13).^(0.5)*exp(betaC3*(-4*DeltaGH/RT+2*F_con*dPsi/RT))*A*B^2-exp((betaC3-1)*(-4*DeltaGH/RT+2*F_con*dPsi/RT))*C*D^2*Keq(13).^(-0.5))/deno;
% %--------------------------------------------------------------
% %----------------------------------------------
% %R 14Complex IV     2CytCred+0.5O2+4Hm---2CytCoxi+H2O+2Hi
% %CytCred+0.25O2+2Hm---CytCoxi+0.5H2O+Hi
betaC4=0.5;
KA14=680e-6;
KB14=5.4e-6;
KC14=79.2e-6;
% KA14=734e-6;
% KB14=1.51e-6;
% KC14=734e-6;
A=CytCred;
B=O2;
C=CytCoxi;
Keq(14)=Keq0(14)*10^(7-pH_m)*10^(7-pH_m);
%-----------------------------------
%  2CytCred+0.5O2+4Hm---2CytCoxi+H2O+2Hi
deno=(1+A^2/KA14^2+C^2/KC14^2)*(1+B^0.5/KB14^0.5);
J_CIV =Vmaxf(14)/KA14^2/KB14^0.5*(Keq(14).^0.5*exp(betaC4*(-2*DeltaGH/RT-2*F_con*dPsi/RT))*A^2*B^0.5-exp((betaC4-1)*(-2*DeltaGH/RT-2*F_con*dPsi/RT))*C^2*Keq(14).^(-0.5))/deno;

%----------------------------------------------
% %ADPm+Pi+3Hi+Hm->ATPm+3Hm
betaC5=0.5;
KA15=100*0.01e-3;% no data
KB15=3e-3;
KC15=100*0.01e-3;
A=ADPm;
B=Pim;
C=ATPm;
Keq(15)=Keq0(15)*10^(7-pH_m);
deno=(1+A/KA15+C/KC15)*(1+B/KB15);
J_CV =Vmaxf(15)/KA15/KB15*(Keq(15).^0.5*exp(3*betaC5*DeltaGH/RT)*A*B-exp(3*(betaC5-1)*DeltaGH/RT)*C*Keq(15).^(-0.5))/deno;
% % %----------------
% Tmax

%% Free Mg and free and Mg-bound ATP and ADP concentrations


%% Transport fluxess
%References: [1] F.Palmieri, Minireview, Mitochondrial carrier proteins
%1994
%-----------------------END
%SUC--Pi   SUCe+Pim -- SUCm+Pie   Antipoter   A1+B2=A2+B1
% K_DCC_MAL=;   ref value =0.7e-3
KDCC_PI=0.93e-3;  %[1]
KDCC_SUC=1e-3; 
% KDCC_MAL=6e-3;
KDCC_MAL=1.17e-3; %average value from ref [1]
%-------------------------------
deno=1+SUCe/KDCC_SUC+SUCm/KDCC_SUC+Pie/KDCC_PI+Pim/KDCC_PI...
    +MALe/KDCC_MAL+MALm/KDCC_MAL...
    +MALe*Pim/KDCC_MAL/KDCC_PI+MALm*Pie/KDCC_MAL/KDCC_PI...
    +SUCe*Pim/KDCC_SUC/KDCC_PI+SUCm*Pie/KDCC_SUC/KDCC_PI;
T_SUC_Pi=2*Tmax(1)/KDCC_SUC/KDCC_PI*(SUCe*Pim-Pie*SUCm)/deno;
%---------------------------------------------------------------
%Pi-MAL   Pim+MALe---Pie+MALme  Antipoter 
deno=1+MALe/KDCC_MAL+MALm/KDCC_MAL+Pie/KDCC_PI+Pim/KDCC_PI...
    +SUCe/KDCC_SUC+SUCm/KDCC_SUC...
    +SUCe*Pim/KDCC_SUC/KDCC_PI+SUCm*Pie/KDCC_SUC/KDCC_PI...
    +MALe*Pim/KDCC_MAL/KDCC_PI+MALm*Pie/KDCC_MAL/KDCC_PI;
T_MAL_Pi=1*Tmax(2)/KDCC_MAL/KDCC_PI*(MALe*Pim-Pie*MALm)/deno;
%-------------------------------------------------------
%aKG-MAL   aKGm +MALe-aKGe+MALm  Antipoter 
 K_aKG=0.24e-3;    %[] ref values=0.31e-3 0.17e-3
 K_MAL=1e-3;   %[]ref values=1.36e-3 0.71e-3
deno=1+MALe/K_MAL+MALm/K_MAL+aKGm/K_aKG+aKGe/K_aKG...
    +aKGm*MALe/K_aKG/K_MAL+aKGe*MALm/K_aKG/K_MAL;
T_MAL_aKG=Tmax(3)/K_aKG/K_MAL*(aKGm*MALe-aKGe*MALm)/deno;
%-----------------------------------------------------------------
%%MALc-HCIT  Antipoter   MALe+Hm+CITm-->He+CITe+MALm  Tricarboxylate 
%%Carrier (TCC)
KTCC_MAL=0.25e-3;   %M(5) REF VALUES:  0.25 and 0.06
KTCC_CIT=1e-3;   %M(6)  16.35e-3
KTCC_H=1e-7;  
deno=1+MALe/KTCC_MAL+MALm/KTCC_MAL+CITm*Hm/KTCC_CIT/KTCC_H+CITe*He/KTCC_CIT/KTCC_H...
    +Hm*CITm*MALe/KTCC_CIT/KTCC_MAL/KTCC_H+He*CITe*MALm/KTCC_CIT/KTCC_MAL/KTCC_H;
T_MAL_HCIT=Tmax(4)/KTCC_MAL/KTCC_CIT/KTCC_H*(MALe*Hm*CITm-MALm*CITe*He)/deno;



%-------------------------------------
%PYR-H co-transporter between e and m (PYRH)  
KPYRH_PYR=0.24e-3;    % M(7)
KPYRH_H=1e-7; 
deno=1+PYRm/KPYRH_PYR+PYRe/KPYRH_PYR+Hm/KPYRH_H+He/KPYRH_H+PYRm*Hm/KPYRH_H/KPYRH_PYR+PYRe*He/KPYRH_H/KPYRH_PYR;
T_PYRH=Tmax(5)/KPYRH_PYR/KPYRH_H*(PYRe*He-PYRm*Hm)/deno;
%----------------------------------------------------------
%Phosphate-H cotransporter (PIC)

KPIC_PI=9.4e-3;
KPIC_H=1e-7;
deno=1+Pim/KPIC_PI+Pie/KPIC_PI+Hm/KPIC_H+He/KPIC_H+Pim*Hm/KPIC_PI/KPIC_H+Pie*He/KPIC_H/KPIC_PI;
T_PIC=Tmax(6)/KPIC_H/KPIC_PI*(Pie*He-Pim*Hm)/deno;
%-------------------------------------------------------------------------
K_ATP=10e-6;
K_ADP=10e-6;
beta_ANT=0.6;
deno=(K_ATP*K_ADP+ADPe+ATPe*exp((beta_ANT-1)*F_con*dPsi/RT))*(K_ATP*K_ADP+ADPm+ATPm*exp(beta_ANT*F_con*dPsi/RT));
T_ANT =Tmax(7)*(exp(beta_ANT*dPsi*F_con/RT)*ADPe*ATPm-exp((beta_ANT-1)*dPsi*F_con/RT)*ADPm*ATPe)/deno;

%---------------------------------------------------------------------
%GLUe*He-GLUm*Hm
KGLU_GLU=1e-3; %?????
KGLU_H=1e-7;
deno=1+GLUe/KGLU_GLU+GLUm/KGLU_GLU+He/KGLU_H+Hm/KGLU_H+GLUe*He/KGLU_H/KGLU_GLU+GLUm*Hm/KGLU_H/KGLU_GLU;
T_GLUH=Tmax(8)/KGLU_H/KGLU_GLU*(GLUe*He-GLUm*Hm)/deno;
%------------------------------------------------
KAG_GLU=0.25e-3;       
KAG_ASP=0.12e-3;
KAG_H=1e-7;
%ASPe +HGLUm----ASPm +HGLUe
deno=1+Hm*GLUm/KAG_GLU/KAG_H+He*GLUe/KAG_GLU/KAG_H+ASPm/KAG_ASP+ASPe/KAG_ASP...
    +Hm*GLUm*ASPe/KAG_GLU/KAG_H/KAG_ASP+He*GLUe*ASPm/KAG_GLU/KAG_H/KAG_ASP;
T_ASP_GLU=Tmax(9)/KAG_GLU/KAG_ASP/KAG_H*(exp(dPsi.*F_con./(2*RT))*GLUm*Hm*ASPe-exp(-dPsi.*F_con./(2*RT))*ASPm*He*GLUe)/deno;
K_Hleak=1e-7;
if dPsi>1e-9
 
     T_HLEAK=Tmax(10)/K_Hleak*(He*exp(dPsi.*F_con./(2*RT))-Hm*exp(-dPsi.*F_con./(2*RT)));

else    %no membrane potential
    T_HLEAK=0;
end
T_FUM_Pi=0;
PS2=0.5*6.88*Vm*60;%Huang et al  6.88+-0.39 moles/(liter mitochondria volume*s)
alfa=F_con/RT;
T_R123=alfa*PS2*dPsi/(exp(alfa*dPsi)-1)*(exp(alfa*dPsi)*R123e-R123m);
%% Metabolic reaction and transport fluxes

%----
rtfluxes=[
J_PDH;
J_CITS;
J_CITD;
J_AKGD;
J_SCAS;
J_NDK;
J_SUCD;
J_FUM;
J_MALD;
J_GOT;
J_CI;
J_CII;
J_CIII;
J_CIV;
J_CV;
T_SUC_Pi;
T_MAL_Pi;
T_MAL_aKG;
T_MAL_HCIT;
T_PYRH;
T_PIC;
T_ANT;
T_GLUH;
T_ASP_GLU;
T_HLEAK;
T_FUM_Pi;
T_R123;];