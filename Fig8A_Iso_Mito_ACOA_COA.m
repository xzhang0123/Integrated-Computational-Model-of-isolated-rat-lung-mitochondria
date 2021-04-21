
clear all
close all
clc
format long
%-----------------------
 size_index=dlmread('text_size.txt');
text_size=size_index(1);
text_size2=size_index(2);
line_width=size_index(3);
marker_size=size_index(4);

%-------------------------
%%  Parameter Setup
%substrates=input('Substrates used in simulation:\n Enter 1 if no substrate;\n Enter 2 for PYR+MAL;\n Enter 3 for SUC;\n Enter 4 for MAL+GLU;\n')
substrates=6;
bufferpH=7.4;
global  Tem C F_con R_con Ve Vm Vi closed_system ROTi AAi...
    iPie  iADPe iATPe iPYRe iMALe iCITe iaKGe iSUCe iFUMe iGLUe iASPe...
    iGLUm iASPm iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iNADm iNADHm...
    iUQm iUQH2m iCytCoxi iCytCred iADPm iATPm iGDPm iGTPm iCOAm iACOAm iPim...
    iFADm iFADH2m iHm iHe idPsi iO2 iR123e iR123m iReBe iRmBm
iPie=1;  iADPe=2; iATPe=3; iPYRe=4; iMALe=5; iCITe=6; iaKGe=7; iSUCe=8; iFUMe=9; iGLUe=10;
iASPe=11; iGLUm=12; iASPm=13; iPYRm=14; iOXAm=15; iCITm=16; iaKGm=17; iSCAm=18; iSUCm=19; iFUMm=20;
iMALm=21; iNADm=22; iNADHm=23; iUQm=24; iUQH2m=25; iCytCoxi=26; iCytCred=27; iADPm=28; iATPm=29;iGDPm=30;
iGTPm=31; iCOAm=32; iACOAm=33; iPim=34; iFADm=35; iFADH2m=36; iHm=37; iHe=38; idPsi=39; iO2=40;
%----Test:include R123 as state variables
iR123e=41; iR123m=42; iReBe=43; iRmBm=44;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
 Tem=303.15; %K      30 oC

R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
ROTi=1;
AAi=1;
%volumes
Ve  =   2/5*1e-3; %L 
Vm  =1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)%---------------------
closed_system=0;% open system, oxygen is constant
%---------------------------------



Q10=2.5*ones(25,1);
Tem_standard=303.15;  %K  30
Q10_factor=Q10.^((Tem-Tem_standard)/10);

Para=ones(27,1);

 Para(1:25)=1*Q10_factor.*Para(1:25);
  %%  Define t_step and t_final
t_step      =   0.1;   %min
%% Run Simulation
   Ti0=1; %time before measuring
time1=10; %Length of State 2  (time before adding ADP)
tic
MALstep=[0:1e-3:10e-3];
Pistep=[0:0.5e-3:4e-3];
%---------------
CIT_rel1=zeros(time1/t_step+1,5);
PYR_upt1=zeros(time1/t_step+1,5);
MAL_upt1=zeros(time1/t_step+1,5);

%-------------
for i=1:1:length(MALstep)
    IC=Set_Initial_Concentrations(substrates,bufferpH);
   IC(iMALe)=MALstep(i);
 IC(iPYRe)=5e-3;
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on');
%       options = odeset('NonNegative',[1:42]);
[T00,C00] = ode15s(@odeq,[0:t_step:Ti0],IC,options,substrates,Para,2);
IC1=C00(end,:);
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,substrates,Para,2);
C=C1;

Tfluxes=zeros(length(C(:,1)),12)';
Rfluxes=zeros(length(C(:,1)),15)';
for istep=1:1:(length(C(1:end,1)))
    RTfluxes(:,istep)=fluxes(C(istep,:),Para);
end

Rfluxes=RTfluxes(:,1:15);
Tfluxes=RTfluxes(:,16:27);
% PYRm1(:,i)=1e3*C(:,iPYRm);
ACOAm1=1e3*C(end,iACOAm);
 COAm1=1e3*C(end,iCOAm);
 C_ACOA(i)=ACOAm1;
 C_COA(i)=COAm1;
ratio(i)=ACOAm1/COAm1;
% clear T1 C1 C T00 C00 T
end

   MAL_index=[0 5e-3];
   COA=[0.32 0.69];
   SE_COA=0.15*COA;
   ACOA=[1.37 0.53];
   SE_ACOA=0.15*ACOA;
 mean_COA_ratio=[0.45/0.74 0.32/1.37 0.53/0.69];

set(figure(3),'Units','inches','Position',[2 1 5 4])

figure(3)
h1=plot(1e3*MALstep,C_ACOA./(C_ACOA+C_COA),'K',1e3*MALstep,C_COA./(C_ACOA+C_COA),'K','linewidth',2)
xlabel('MAL (mM)')
ylabel('%ACOA and %COA')
text(6,0.5,'[PYR] = 5mM','Fontsize',text_size2)
xlim([-0.5 10])
box off
title('%ACOA and %COA')
 set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
hold on
h2=errorbar(1e3*MAL_index,COA./(ACOA+COA),SE_COA./(ACOA+COA),'^','MarkerSize',12,'Color','b','linewidth',2,'CapSize',10)
hold on
h3=errorbar(1e3*MAL_index,ACOA./(ACOA+COA),SE_ACOA./(ACOA+COA),'o','MarkerSize',12,'Color','r','linewidth',2,'CapSize',10)
hold off
legend([h2 h3 h1(2) ],'COA','ACOA','Model')
legend boxoff