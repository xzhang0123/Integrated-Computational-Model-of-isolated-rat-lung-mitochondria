
function  IC=Set_Initial_Concentrations(substrates,bufferpH)
%%  Set initial conditions (concentration, [M] unless otherwise noted)
%% mM=10^-6 mol/ml =10^-3 mol/dm3
%10.^-3 mol/dm3 =10.^-6 mol/cm3
global k0 Tem F_con R_con Ve Vm Vi ROTi AAi closed_system...
 iPie  iADPe iATPe iPYRe iMALe iCITe iaKGe iSUCe iFUMe iGLUe iASPe...
    iGLUm iASPm iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iNADm iNADHm...
    iUQm iUQH2m iCytCoxi iCytCred iADPm iATPm iGDPm iGTPm iCOAm iACOAm iPim...
    iFADm iFADH2m iHm iHe idPsi iO2 iR123e iR123m iReBe iRmBm

format long

%NADH_M=0.87e-3; % total NADH+NAD
% NADH_M=3e-3; % total NADH+NAD
% Total_CytC=0.31e-3;
Total_TCA_Pool=7e-3;
Total_NADH=1.73e-9/Vm;  %mass/Volume,  UNIT: mol/L
Total_FAD= 0.7e-3;    %ref= Sonia Cortassa et al. Mitochondria -- computational study.
% Total_CytC=0.33e-9/(Vm+Vi)*0.65;%
Total_CytC=0.33e-9/Vi;
Total_ANP=6.4e-9/Vm;
Total_COA=0.925e-9/Vm;
UQpool=0.52e-3;
Pie0 = 4e-3;        %1
ADPe0=0;           %2
ATPe0 = 0.000;       %3
PYRe0 = 0;        %4

MALe0=0e-3;     %5
CITe0 =0;  %6
aKGe0 =0;    %7
SUCe0 = 0;    %8
FUMe0=0;
GLUe0 = 0.0;     %9
ASPe0 =0e-3;      %10
GLUm0 =0e-3;     %11
ASPm0 =12e-3;       %12
PYRm0 =3e-3;        %13     
%OAA0 = 0.001e-3;      %14
OAA0=Total_TCA_Pool/7;
CITm0 =Total_TCA_Pool/7;        %15
aKGm0 = Total_TCA_Pool/7;        %16
COApool=Total_COA;
SCA0 = Total_TCA_Pool/7;       %17
SUCm0 = Total_TCA_Pool/7;        %18
FUMm0=Total_TCA_Pool/7;
MALm0 =Total_TCA_Pool/7;        %19
NADm0=0.25*Total_NADH; %20
NADHm0=0.75*Total_NADH;  %21

UQm0=0.5*UQpool;   %22
UQH2m0=0.5*UQpool;      %23
% CytCoxi0 = 0.9*2.7e-3;  %24
% CytCred0 = 0.1*2.7e-3;    %25
 CytCoxi0 = 0.9*Total_CytC;  %24
 CytCred0 = 0.1*Total_CytC;  %25
% ADPm0 = 0.2e-3;   %26
% ATPm0 =8e-3;       %27
ADPm0 = 0.025*Total_ANP;   %26
ATPm0 =0.975*Total_ANP;       %27
GDPm0 =0.1*ADPm0;   %28
GTPm0 =0.1*ATPm0;       %29
COA0 = 0.5*(COApool);       %30
ACOA0 = 0.5*(COApool);     %31
Pim0 = 4e-3;          %32
FAD0 =0.5* Total_FAD;        %33  FAD content is assumed to be the same as NADH content
FADH20 = 0.5*Total_FAD;      %34
Hm0 = 10^(-7.6);    %35
He0 = 10^(-bufferpH); %36
deltaFi0 = 150;      %37 millivolts
% Membrane potential between intermembrane space and mito
O20 = 0.296e-3;        % 40
R123e0=0e-9;       %39
R123m0=0;
ReBe0=0;    %40
RmBm0=0;



%% choose substrates
if substrates==1  %no substrate
  
elseif substrates==2  %PYR+MAL
 CITe0=5e-3;    %14
 PYRe0=5e-3;
elseif substrates==3  %SUC
SUCe0=5e-3;

elseif substrates==4  %both
SUCe0=5e-3;
MALe0 =5e-3; 
PYRe0 =5e-3; 
elseif substrates==5 %5mM glu+mal
  
  MALe0 =5e-3; 
  GLUe0 = 5e-3;   
elseif substrates==6 %5mM PYR and 5mM MAL 
    
MALe0 =5e-3; 
PYRe0 =5e-3; 

elseif substrates==7 %
    %-------------
    CITe0=5e-3;    %14
 %-----------------------
elseif substrates==8 %fully oxidized 
    
MALe0 =5e-3; 
PYRe0 =5e-3; 
Pim0=0;
NADm0=1*NADH_M; %20
NADHm0=0*NADH_M;  %21
UQpool=1.35e-3;
UQm0=1*UQpool;   %22
UQH2m0=0*UQpool;      %23
CytCoxi0 = 1*2.7e-3;  %24
CytCred0 = 0*2.7e-3;    %25
FAD0 =1* 0.122e-3;        %33
FADH20 = 0*0.122e-3;      %34
%-----


%  aKGm0=5e-3;
%  MALm0=7e-3;
%  SUCm0=5e-3;
elseif substrates==9 %CIT+PYR
    
CITe0 =5e-3; 
PYRe0 =5e-3; 
elseif substrates==10 %aKG+PYR
    
aKGe0 =5e-3; 
PYRe0 =5e-3; 
elseif substrates==11 %PYR
     
PYRe0 =5e-3; 
elseif substrates==12 %MAL
     
MALe0 =5e-3; 
elseif substrates==13 %GLU
     
GLUe0 =10e-3; 
elseif substrates==14 %akg
     
aKGe0 =3e-3; 
end
% emf
%  DeltaGH=Faraday_constant*deltaFi0+ R_constant *Tem*log(Hi0/Hm0);
%% save initialconcentrations
 initialconditions = [

    Pie0; %1
   
 
    ADPe0;  %2
    ATPe0;  %3
 PYRe0;%4
    MALe0;  %5
    CITe0;  %6
    aKGe0;   %7
    SUCe0;   %8
    FUMe0;
    GLUe0;   %10
    ASPe0     %11
    GLUm0;    %2
    ASPm0;    %13
    PYRm0;  %14
    OAA0;    %15
    CITm0;    %16
    aKGm0;   %17
    SCA0;  %18
    SUCm0;  %19
    FUMm0;
    MALm0;  %21
    NADm0;  %22
    NADHm0; %23)
    UQm0;   %24
    UQH2m0;     %25
    CytCoxi0;    %26
    CytCred0;    %27
    ADPm0;   %28
    ATPm0;  %29
    GDPm0; %30
    GTPm0; %31
    COA0;   %32
    ACOA0;  %33
    Pim0;   %34
    FAD0;   %35
    FADH20;%36
    Hm0;%37
    He0;%38
    deltaFi0;%39
    O20;  %40
     R123e0;%41
     R123m0; %42 
     ReBe0;
    RmBm0;

   ];%44
IC=initialconditions;
