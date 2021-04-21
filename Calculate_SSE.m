function SSE=Calculate_SSE(Para,Edata)
format long
%Objective function for TCA cycle data at state 2
%%  Parameter Setup
global  Tem F_con R_con Ve Vm Vi ROTi AAi closed_system...
 iPie  iADPe iATPe iPYRe iMALe iCITe iaKGe iSUCe iFUMe iGLUe iASPe...
    iGLUm iASPm iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iNADm iNADHm...
    iUQm iUQH2m iCytCoxi iCytCred iADPm iATPm iGDPm iGTPm iCOAm iACOAm iPim...
    iFADm iFADH2m iHm iHe idPsi iO2 iR123e iR123m iReBe iRmBm

%% Run Simulation 1 TCA cycle data
closed_system=0;
substrates=6;
Pi_buffer=5e-3;
Ve  =   1/0.7*1e-3; %L
bufferpH=7.4;
% Tem=303.15; %K      30 oC
PYRstep=[0 0.5e-3 1e-3 5e-3 10e-3];
MALstep=[0 0.5e-3 1e-3 5e-3 10e-3];
for i=1:1:5
    IC=Set_Initial_Concentrations(substrates,bufferpH);
   IC(iPYRe)=PYRstep(i);
   
time1=1;
 time2=10; 
 %options = odeset('MaxStep',0.1,'NonNegative',[1:44]);
 options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
           'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:44]);
[~,C00] = ode15s(@odeq,[0:0.2:time1],IC,options,substrates,Para,2);
IC1_dPYR=C00(end,:);

[~,C_dPYR] = ode15s(@odeq,[0:2:10],IC1_dPYR,options,substrates,Para,2);
CIT_rel1(:,i)=(C_dPYR(:,iCITe)-C_dPYR(1,iCITe))*Ve*1e9; %nmol
PYR_upt1(:,i)=(C_dPYR(1,iPYRe)-C_dPYR(:,iPYRe))*Ve*1e9;
MAL_upt1(:,i)=(C_dPYR(1,iMALe)-C_dPYR(:,iMALe))*Ve*1e9;
end
CIT1=[CIT_rel1(2:end,1);CIT_rel1(2:end,2);CIT_rel1(2:end,3);CIT_rel1(2:end,4);CIT_rel1(2:end,5)];
PYR1=[PYR_upt1(2:end,1);PYR_upt1(2:end,2);PYR_upt1(2:end,3);PYR_upt1(2:end,4);PYR_upt1(2:end,5)]; %PYR=0.5 and 10
MAL1=[MAL_upt1(2:end,1);MAL_upt1(2:end,2);MAL_upt1(2:end,3);MAL_upt1(2:end,4);MAL_upt1(2:end,5)];

for i=1:1:5
    IC=Set_Initial_Concentrations(substrates,bufferpH);
   IC(iMALe)=MALstep(i);
[~,C00] = ode15s(@odeq,[0:0.2:time1],IC,options,substrates,Para,2);
IC1_dMAL=C00(end,:);
[~,C_dMAL] = ode15s(@odeq,[0:2:time2],IC1_dMAL,options,substrates,Para,2);

CIT_rel2(:,i)=(C_dMAL(:,iCITe)-C_dMAL(1,iCITe))*Ve*1e9; %nmol
PYR_upt2(:,i)=(C_dMAL(1,iPYRe)-C_dMAL(:,iPYRe))*Ve*1e9;
MAL_upt2(:,i)=(C_dMAL(1,iMALe)-C_dMAL(:,iMALe))*Ve*1e9;
end
CIT2=[CIT_rel2(2:end,1);CIT_rel2(2:end,2);CIT_rel2(2:end,3);CIT_rel2(2:end,4);CIT_rel2(2:end,5)];
PYR2=[PYR_upt2(2:end,1);PYR_upt2(2:end,2);PYR_upt2(2:end,3);PYR_upt2(2:end,4);PYR_upt2(2:end,5)];
MAL2=[MAL_upt2(2:end,1);MAL_upt2(2:end,2);MAL_upt2(2:end,3);MAL_upt2(2:end,4);MAL_upt2(2:end,5)];

Sim_TCAdata1=[CIT1' PYR1' MAL1' CIT2' PYR2' MAL2']';
%Simulation 2-----------OCR 4 Susbstrates------
closed_system=1;
t_step=0.5;
Ve  =   1.9/0.4*1e-3; %L  
ADP_amount=300e-6;
substrates=1;
time1=1;
time2=1;
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iCITe)=5e-3;
IC(iPYRe)=5e-3;
IC(iPie)=4e-3;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+1*ADP_amount; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
T=[T1; T2(2:end)+time1;];
C=[C1; C2(2:end,:)];
OC_CITPYR=1e6*C(:,iO2);
%-----------------------------------
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iSUCe)=5e-3;
IC(iPYRe)=5e-3;
IC(iPie)=4e-3;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+1*ADP_amount; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
T=[T1; T2(2:end)+time1;];
C=[C1; C2(2:end,:)];
OC_SUCPYR=1e6*C(:,iO2);
%-----------------------------------
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iaKGe)=5e-3;
IC(iPYRe)=5e-3;
IC(iPie)=4e-3;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+1*ADP_amount; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
T=[T1; T2(2:end)+time1;];
C=[C1; C2(2:end,:)];
OC_AKGPYR=1e6*C(:,iO2);
%-----------------------------------
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iMALe)=5e-3;
IC(iPYRe)=5e-3;
IC(iPie)=4e-3;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+1*ADP_amount; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
T=[T1; T2(2:end)+time1;];
C=[C1; C2(2:end,:)];
OC_MALPYR=1e6*C(:,iO2);
Sim_OCRdata2=[OC_CITPYR; OC_SUCPYR; OC_AKGPYR; OC_MALPYR;];
%----------------Simulation 3--Experiment condition---23degrees OCR------------------------
Ve  =  1*1e-3; %L
bufferpH=7.2;
Pi_buffer=5e-3;
%----Temperature correction
Q10=2.5*ones(25,1);
Tem_standard=303.15;  %K  30
 %Tem_sim2=303.15-7; %K      23 oC
 Tem=303.15-7;
Q10_factor=Q10.^((Tem-Tem_standard)/10);
%------------------------------------------
Para(1:25)=1*Q10_factor.*Para(1:25);
substrates=0; %PYR MAL
time1=3.5;
time2=3;
% time3=1.9;
t_step=0.2;
options = odeset('NonNegative',[1:44]);
IC=Set_Initial_Concentrations(substrates,bufferpH);
Para_PM=Para;
Para_PM(25)=0.40*Para(25);
IC(iPYRe)=10e-3;
IC(iMALe)=5e-3;
IC(iPie)=Pi_buffer;
IC(iO2)=295.4*1e-6;
[~,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para_PM,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)=IC2(iADPe)+ 100e-6; %add ADP, Unit(Molar) 
[~,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para_PM,2);
C=[C1; C2(2:end,:)];
Sim_Oxy_PM=1e6*C(:,iO2);

clear C 
%----------------
IC=Set_Initial_Concentrations(substrates,bufferpH);
Para_SUC=Para;
Para_SUC(25)=0.40*Para(25);
IC(iSUCe)=7e-3;
IC(iPie)=Pi_buffer;
IC(iO2)=295.4*1e-6;
[~,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para_SUC,2);

IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)=IC2(iADPe)+ 100e-6; %add ADP, Unit(Molar) 
[~,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para_SUC,2);
% 
% IC3=C2(end,:);
%  Para_SUC(25)= UP_factor*Para_SUC(25); %leak
% [~,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para_SUC,2);
% C=[C1; C2(2:end,:);C3(2:end,:)];

C=[C1; C2(2:end,:)];
Sim_Oxy_SUC=1e6*C(:,iO2);
Sim_data=[Sim_TCAdata1./length(Sim_TCAdata1);Sim_Oxy_PM/length(Sim_Oxy_PM);Sim_Oxy_SUC./length(Sim_Oxy_SUC);Sim_OCRdata2./length(Sim_OCRdata2)];
% length(Sim_data)
% length(Sim_TCAdata1)
% length(Sim_Oxy_PM)
% length(Sim_Oxy_SUC)
% length(Sim_OCRdata2)
SSE=sum((Sim_data-Edata).^2)
% size(Sim_ydata)
% size(Sim_Oxy_SUC)