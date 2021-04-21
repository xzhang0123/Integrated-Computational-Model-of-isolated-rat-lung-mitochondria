
clear all
close all
clc
format long
%-----------------------
%figure 1&2 experimental data and fits,%2 by 3
%Figure 3 time course of species (different PYR concentrations)


%experiment conditions :
%pH=7.4   Pi=4mM
%mitchondria =0.7mg,temperature =30oC
text_size=16;
mk_size=12;
line_width=1.5;
%-------------------------
%%  Parameter Setup
%substrates=input('Substrates used in simulation:\n Enter 1 if no substrate;\n Enter 2 for PYR+MAL;\n Enter 3 for SUC;\n Enter 4 for MAL+GLU;\n')
substrates=1;
bufferpH=7.2;
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
 %Tem=303.15; %K      30 oC
  Tem=298.15; %K      25 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
ROTi=1;
AAi=1;
%volumes
Ve  =   1/2.5*1e-3; %L
Vm  = 1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
%---------------------
closed_system=0;% 0: open system, oxygen is constant; 1: closed system
%---------------------------------
IC=Set_Initial_Concentrations(substrates,bufferpH);

 IC(iPie)=5e-3;
% Temperature correction---------------
Q10=2.5*ones(25,1);
Tem_standard=303.15;  %K  30
 Q10_factor=Q10.^((Tem-Tem_standard)/10);

Para=ones(27,1);
 Para(1:25)=1*Q10_factor.*Para(1:25);
Para(25)=1*Para(25); %LEAK
  %%  Define t_step and t_final
t_step      =   0.1;   %min
%% Run Simulation
   Ti0=1; %time before measuring

tic


    time0=2;
    time1=3;
    time2=1.2;
    time3=13;
 IC=Set_Initial_Concentrations(substrates,bufferpH);

options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on');
%       options = odeset('NonNegative',[1:42]);
[T00_pre,C00_pre] = ode15s(@odeq,[0:t_step:time0],IC,options,substrates,Para,2);
IC0=C00_pre(end,:);
[T00,C00] = ode15s(@odeq,[0:t_step:Ti0],IC0,options,substrates,Para,2);
IC1=C00(end,:);
IC1(iPYRe)= 2.5e-3; %add PYR, Unit(Molar) 
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,substrates,Para,2);

IC2=C1(end,:); %Intial concentration of State III
IC2(iMALe)= 2.5e-3; %add MAL, Unit(Molar) 
%IC2(iPYRe)= 2.5e-3;
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);

IC3=C2(end,:); %Intial concentration of State III
IC3(iADPe)= 200e-6; %add ADP, Unit(Molar) 
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para,2);

T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2;];
C=[C1; C2(2:end,:);C3(2:end,:)];
Tfluxes=zeros(length(C(:,1)),12)';
Rfluxes=zeros(length(C(:,1)),15)';
for istep=1:1:(length(C(1:end,1)))
    RTfluxes(:,istep)=fluxes(C(istep,:),Para);
end
Rfluxes=RTfluxes(:,1:15);
Tfluxes=RTfluxes(:,16:27);
PYRm1=1e3*C(:,iPYRm);
ACOAm1=1e3*C(:,iACOAm);
COAm1=1e3*C(:,iCOAm);
OXAm1=1e3*C(:,iOXAm);
CITm1=1e3*C(:,iCITm);
aKGm1=1e3*C(:,iaKGm);
SCAm1=1e3*C(:,iSCAm);
SUCm1=1e3*C(:,iSUCm);
FUMm1=1e3*C(:,iFUMm);
MALm1=1e3*C(:,iMALm);
OXAm1=1e3*C(:,iOXAm);
NADHm1=1e3*C(:,iNADHm);
NADm1=1e3*C(:,iNADm);
FADH2m1=1e3*C(:,iFADH2m);
UQm1=1e3*C(:,iUQm);
CytCredm1=1e3*C(:,iCytCred);
dPsi1=C(:,idPsi);
UQH2m1=1e3*C(:,iUQH2m);
ATPm1=1e3*C(:,iATPm);
%------------------
PYRDH1=1e9*Rfluxes(1,:);
CITSYN1=1e9*Rfluxes(2,:);
CITDH1=1e9*Rfluxes(3,:);
AKGDH1=1e9*Rfluxes(4,:);
SCAS1=1e9*Rfluxes(5,:);
SUCDH1=1e9*Rfluxes(7,:);
FUM1=1e9*Rfluxes(8,:);
MALDH1=1e9*Rfluxes(9,:);
CI1=1e9*Rfluxes(11,:);
CII1=1e9*Rfluxes(12,:);
CIII1=1e9*Rfluxes(13,:);
CIV1=1e9*Rfluxes(14,:);
CV1=1e9*Rfluxes(15,:);
DCC1=1e9*Tfluxes(1,:);
DCC2=1e9*Tfluxes(2,:);
AKGMAL1=1e9*Tfluxes(3,:);
TCC1=1e9*Tfluxes(4,:);
PYRH1=1e9*Tfluxes(5,:);
PIC1=1e9*Tfluxes(6,:);
LEAK1=1e9*Tfluxes(10,:);
%-----

%-----------------------
clear T1 C1 C T00 C00 

toc
clear t t_step t_final

  NADH_data=dlmread('NADH_data.txt');
        time_NADH=NADH_data(:,1);
        I_NADH=NADH_data(:,2);
        I_NADH=I_NADH/max(I_NADH)/1.5;
        figure(2)
        plot(time_NADH,I_NADH,'*')
        title('NADH(Experiment)')
        ylim([0 1])
%--------------------------

set(figure(1),'Units','inches','Position',[2 1 5 4])
figure(1)
plot(time_NADH,I_NADH,'b*',T-2,NADHm1./(NADHm1+NADm1),'k','LineWidth',line_width,'MarkerSize',10)
legend('Data','Model')
legend boxoff
xlabel('Time (min)')
ylabel('Normalized Intensity')
 xlim([0.5 4.2])
  ylim([0 1])
  box off 
  set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 title('NADH ratio')
% figure(2)
% plot(time_NADH,I_NADH,'b*','LineWidth',line_width,'MarkerSize',8)
% legend('Model','Data')
% legend boxoff
% xlabel('Time (min)')
% ylabel('Normalized Intensity')
%  xlim([0 4.2])
%   ylim([0 1])
%   box off 
%   set(gcf,'color','w')
% set(gca,'Fontsize',text_size,'LineWidth',line_width)
% title('NADH ratio')
% 
%       figure(3)
% plot(T-2,NADHm1./(NADHm1+NADm1),'b','LineWidth',line_width,'MarkerSize',8)
% legend('Model','Data')
% legend boxoff
% xlabel('Time (min)')
% ylabel('Normalized Intensity')
%  xlim([0 4.2])
%   ylim([0 1])
%   box off 
%   set(gcf,'color','w')
% set(gca,'Fontsize',text_size,'LineWidth',line_width)
% title('NADH ratio')
%        