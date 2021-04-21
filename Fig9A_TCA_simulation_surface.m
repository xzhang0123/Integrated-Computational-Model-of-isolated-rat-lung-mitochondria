
clear all
close all
clc
format long
%-----------------------

text_size=16;
mk_size=12;
line_width=1.5;
%experiment conditions :
%pH=7.4   Pi=4mM
%mitchondria =0.7mg,temperature =30oC
%-------------------------
%%  Parameter Setup
%substrates=input('Substrates used in simulation:\n Enter 1 if no substrate;\n Enter 2 for PYR+MAL;\n Enter 3 for SUC;\n Enter 4 for MAL+GLU;\n')
substrates=1;
bufferpH=7.2;
global F_con R_con Tem Ve Vm Vi ROTi AAi closed_system...
    iPie  iADPe iATPe iPYRe iMALe iCITe iaKGe iSUCe iFUMe iGLUe iASPe...
    iGLUm iASPm iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iNADm iNADHm...
    iUQm iUQH2m iCytCoxi iCytCred iADPm iATPm iGDPm iGTPm iCOAm iACOAm iPim...
    iFADm iFADH2m iHm iHe idPsi iO2 iR123e iR123m iReBe iRmBm
%-----------
%define index
iPie=1;  iADPe=2; iATPe=3; iPYRe=4; iMALe=5; iCITe=6; iaKGe=7; iSUCe=8; iFUMe=9; iGLUe=10;
iASPe=11; iGLUm=12; iASPm=13; iPYRm=14; iOXAm=15; iCITm=16; iaKGm=17; iSCAm=18; iSUCm=19; iFUMm=20;
iMALm=21; iNADm=22; iNADHm=23; iUQm=24; iUQH2m=25; iCytCoxi=26; iCytCred=27; iADPm=28; iATPm=29;iGDPm=30;
iGTPm=31; iCOAm=32; iACOAm=33; iPim=34; iFADm=35; iFADH2m=36; iHm=37; iHe=38; idPsi=39; iO2=40;
%----Test:include R123 as state variables
iR123e=41; iR123m=42; iReBe=43; iRmBm=44;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
 Tem=303.15; %K      30 oC
  %Tem=298.15; %K      25 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
ROTi=1;
AAi=1;
%volumes
Ve  =   1/0.7*1e-3; %L
%Ve  =   1e-3; %L
%Vmito=3.67e-6;% total volume  L (for 1mg mitochondria) 
Vm  = 1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
Volumes = [Ve; Vm; Vi];
%---------------------
closed_system=0;% open system, oxygen is constant
%---------------------------------
IC=Set_Initial_Concentrations(substrates,bufferpH);



%--------------------------------
% Temperature correction
Q10=2.5*ones(25,1);
Tem_standard=303.15;  %K  30
Q10_factor=Q10.^((Tem-Tem_standard)/10);
%Para = dlmread('EstPara_final1.txt'); 
Para=ones(27,1);
Para(1:25)=1*Q10_factor.*Para(1:25);

%dlmwrite('EstPara_test.txt', Para,'precision','%10.14d');
  %%  Define t_step and t_final
t_step      =   1;   %min
%% Run Simulation
   Ti0=2;
time1=10; %Length of State 2  (time before adding ADP)
tic
PYRstep=[ 0.1e-3 :2e-3 : 10.1e-3];
MALstep=[ 0.1e-3 :2e-3 : 10.1e-3];

for j=1:1:length(PYRstep)
   IC(iPYRe)=PYRstep(j);
for i=1:1:length(MALstep)
   
   IC(iMALe)=MALstep(i);
   IC(iPie)=4e-3;
% options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
%           'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
%           'BDF','on');
      options = odeset('NonNegative',[1:42]);
[T00,C00] = ode15s(@odeq,[0:t_step:Ti0],IC,options,substrates,Para,2);
IC1=C00(end,:);
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,substrates,Para,2);
C=C1;
T=T1;

Tfluxes=zeros(length(C(:,1)),12)';
Rfluxes=zeros(length(C(:,1)),15)';
for istep=1:1:length(C(:,1))
[RTfluxes(:,istep)]=fluxes(C(istep,:),Para);
end
Rfluxes=RTfluxes(1:15,:);
Tfluxes=RTfluxes(16:27,:);
%ETC fluxes and leak
CI(i,j)=1e9*Rfluxes(11,end);
CII(i,j)=1e9*Rfluxes(12,end);
CIII(i,j)=1e9*Rfluxes(13,end);
CIV(i,j)=1e9*Rfluxes(14,end);
CV(i,j)=1e9*Rfluxes(15,end);
leak(i,j)=1e9*Tfluxes(10,end);
%Transport fluxes
DCCSUC(i,j)=1e9*Tfluxes(1,end);
DCCMAL(i,j)=1e9*Tfluxes(2,end);
MALAKG(i,j)=1e9*Tfluxes(3,end);
TCC(i,j)=1e9*Tfluxes(4,end);
PYRH(i,j)=1e9*Tfluxes(5,end);
%PIC(i,j)=Tfluxes(6,end);
%TCA fluxes
PYRDH(i,j)=1e9*Rfluxes(1,end);
CITSYN(i,j)=1e9*Rfluxes(2,end);
CITD(i,j)=1e9*Rfluxes(3,end);
AKGDH(i,j)=1e9*Rfluxes(4,end);
SCAS(i,j)=1e9*Rfluxes(5,end);
%
SUCDH(i,j)=1e9*Rfluxes(7,end);
FUM(i,j)=1e9*Rfluxes(8,end);
MALDH(i,j)=1e9*Rfluxes(9,end);
%---------------------------------
PYRm(i,j)=1e3*C(end,iPYRm);
ACOAm(i,j)=1e3*C(end,iACOAm);
COAm(i,j)=1e3*C(end,iCOAm);
OXAm(i,j)=1e3*C(end,iOXAm);
CITm(i,j)=1e3*C(end,iCITm);
aKGm(i,j)=1e3*C(end,iaKGm);
SCAm(i,j)=1e3*C(end,iSCAm);
SUCm(i,j)=1e3*C(end,iSUCm);

FUMm(i,j)=1e3*C(end,iFUMm);
MALm(i,j)=1e3*C(end,iMALm);
NADHm(i,j)=1e3*C(end,iNADHm);
NADm(i,j)=1e3*C(end,iNADm);
FADH2m(i,j)=1e3*C(end,iFADH2m);
FADm(i,j)=1e3*C(end,iFADm);
UQH2m(i,j)=1e3*C(end,iUQH2m);
UQm(i,j)=1e3*C(end,iUQm);
CytCredm(i,j)=1e3*C(end,iCytCred);
CytCoxim(i,j)=1e3*C(end,iCytCoxi);
dPsi(i,j)=C(end,idPsi);
ATPm(i,j)=1e3*C(end,iATPm);
ADPm(i,j)=1e3*C(end,iADPm);
%------------
clear T1 C1 C T00 C00 T Rfluxes Tfluxes
%IC=Set_Initial_Concentrations(substrates,bufferpH);
end
end

toc
clear t t_step t_final

%-------------PYR-----------
set(figure(1),'Units','inches','Position',[2 1 16 16]) 
figure(1)
subplot(4,4,1)
mesh(1e3*MALstep,1e3*PYRstep,PYRm)
zlim([0 12])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-------------------------------
subplot(4,4,2)
mesh(1e3*MALstep,1e3*PYRstep,ACOAm)
zlim([0 1])

set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%--------------------------------
% subplot(4,4,3)
% mesh(1e3*MALstep,1e3*PYRstep,COAm)
% set(gcf,'color','w')
% set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,3)
mesh(1e3*MALstep,1e3*PYRstep,OXAm)
% zlim([0 1])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,4)
mesh(1e3*MALstep,1e3*PYRstep,CITm)
% zlim([0 40])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,5)
mesh(1e3*MALstep,1e3*PYRstep,aKGm)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,6)
mesh(1e3*MALstep,1e3*PYRstep,SCAm)
zlim([0 0.01])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,7)
mesh(1e3*MALstep,1e3*PYRstep,SUCm)
zlim([0 0.1])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,8)
mesh(1e3*MALstep,1e3*PYRstep,FUMm)
zlim([0 4])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,9)
mesh(1e3*MALstep,1e3*PYRstep,MALm)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,10)
mesh(1e3*MALstep,1e3*PYRstep,NADHm)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,11)
mesh(1e3*MALstep,1e3*PYRstep,FADH2m)
% zlim([0 0.05])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,12)
mesh(1e3*MALstep,1e3*PYRstep,UQH2m)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,13)
mesh(1e3*MALstep,1e3*PYRstep,CytCredm)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,14)
mesh(1e3*MALstep,1e3*PYRstep,ATPm)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)

%-----------------------------------------------------
subplot(4,4,15)
mesh(1e3*MALstep,1e3*PYRstep,dPsi)
zlim([150 180])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%% 
%----------------------------------------------------------
 set(figure(2),'Units','inches','Position',[2 1 16 16]) 
figure(2)
subplot(4,4,1)
mesh(1e3*MALstep,1e3*PYRstep,PYRH)
zlim([0 15])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-------------------------------
subplot(4,4,2)
mesh(1e3*MALstep,1e3*PYRstep,TCC)
zlim([0 15])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%--------------------------------
subplot(4,4,3)
mesh(1e3*MALstep,1e3*PYRstep,DCCMAL)
%zlim([0 50])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,4)
mesh(1e3*MALstep,1e3*PYRstep,DCCSUC)
zlim([-2 0])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)

subplot(4,4,5)
mesh(1e3*MALstep,1e3*PYRstep,PYRDH)
zlim([0 15])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,6)
mesh(1e3*MALstep,1e3*PYRstep,CITSYN)
zlim([0 15])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,7)
mesh(1e3*MALstep,1e3*PYRstep,CITD)
zlim([0 1])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,8)
mesh(1e3*MALstep,1e3*PYRstep,AKGDH)

set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
% subplot(4,4,9)
% mesh(1e3*MALstep,1e3*PYRstep,SCAS)
% zlim([0 5])
% set(gcf,'color','w')
% set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,9)
mesh(1e3*MALstep,1e3*PYRstep,SUCDH)
%zlim([-10 10])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,10)
mesh(1e3*MALstep,1e3*PYRstep,FUM)
%zlim([-10 10])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,11)
mesh(1e3*MALstep,1e3*PYRstep,MALDH)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,12)
mesh(1e3*MALstep,1e3*PYRstep,CI)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,13)
mesh(1e3*MALstep,1e3*PYRstep,CIII)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,14)
mesh(1e3*MALstep,1e3*PYRstep,CIV)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------
subplot(4,4,15)
mesh(1e3*MALstep,1e3*PYRstep,leak)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
%-----------------------------------------------------

