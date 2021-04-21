
clear all
clc
format long

%----
text_size=15;
mk_size=8;
line_width=1.5;
%----------------------
%%  Parameter Setup
bufferpH=7.2;
global  Tem F_con R_con Ve Vm Vi ROTi AAi closed_system...
 iPie  iADPe iATPe iPYRe iMALe iCITe iaKGe iSUCe iFUMe iGLUe iASPe...
    iGLUm iASPm iPYRm iOXAm iCITm iaKGm iSCAm iSUCm iFUMm iMALm iNADm iNADHm...
    iUQm iUQH2m iCytCoxi iCytCred iADPm iATPm iGDPm iGTPm iCOAm iACOAm iPim...
    iFADm iFADH2m iHm iHe idPsi iO2 iR123e iR123m iReBe iRmBm
iPie=1;  iADPe=2; iATPe=3; iPYRe=4; iMALe=5; iCITe=6; iaKGe=7; iSUCe=8; iFUMe=9; iGLUe=10;
iASPe=11; iGLUm=12; iASPm=13; iPYRm=14; iOXAm=15; iCITm=16; iaKGm=17; iSCAm=18; iSUCm=19; iFUMm=20;
iMALm=21; iNADm=22; iNADHm=23; iUQm=24; iUQH2m=25; iCytCoxi=26; iCytCred=27; iADPm=28; iATPm=29;iGDPm=30;
iGTPm=31; iCOAm=32; iACOAm=33; iPim=34; iFADm=35; iFADH2m=36; iHm=37; iHe=38; idPsi=39; iO2=40;
%Test:include R123 as state variables
iR123e=41; iR123m=42; iReBe=43; iRmBm=44;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
  Tem=298.15-2; %K      23 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 Ve  =  1*1e-3; %L    
Vm  = 1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
ROTi=1;
AAi=1;
Volumes = [Ve; Vm; Vi];
closed_system=0;% 0=open system
% % Temperature correction
Q10=2.5*ones(25,1);
Tem_standard=303.15;  %K  30
Q10_factor=Q10.^((Tem-Tem_standard)/10);
Para=ones(27,1);

 Para(1:25)=1*Q10_factor.*Para(1:25);
Para(25)=0.40*Para(25); %Membrane leak in our experiments is ~40% of Evans'

  %%  Define t_step and t_final
t_step      =   0.01;   %min
%% Run Simulation
time0=5; %equalibrating time
time1=3.8; %Length of State 2  (time before adding ADP)
time2=8.1;
time3=8;
substrates=6;
ADP_add=100e-6;
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iMALe)=5e-3;
IC(iPYRe)=10e-3;
IC(iR123e)=0;
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:44]);
[T0,C0] = ode15s(@odeq,[0:t_step:time0],IC,options,substrates,Para,2);
IC1=C0(end,:);
IC1(iR123e)=200e-9; %add dye
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= ADP_add; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
IC3=C2(end,:);
IC3(iADPe)=0.5*ADP_add; %add ADP, Unit(Molar) 
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para,3);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C_PM=[C1; C2(2:end,:);C3(2:end,:)];
Dfi1=C_PM(:,idPsi);
NADH1=C_PM(:,iNADHm);
UQH2m1=C_PM(:,iUQH2m);
CytoCm1=C_PM(:,iCytCred);
Oxygen1=C_PM(:,iO2);
ADPe1=C_PM(:,iADPe);
ADPm1=C_PM(:,iADPm);
R123e_PM=C_PM(:,iR123e);
ReBe_PM=C_PM(:,iReBe);
R123m_PM=C_PM(:,iR123m);
clear C IC1 IC2 IC3
%%----------------------
%SUC
substrates=3;
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iSUCe)=7e-3;
IC(iR123e)=0;
[T0,C0] = ode15s(@odeq,[0:t_step:time0],IC,options,substrates,Para,2);
IC1=C0(end,:);
IC1(iR123e)=200e-9; %add dye
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,substrates,Para,2);

IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= ADP_add; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
IC3=C2(end,:);
IC3(iADPe)= 0.5*ADP_add; %add ADP, Unit(Molar) 
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para,3);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C_SUC=[C1; C2(2:end,:);C3(2:end,:)];
Dfi2=C_SUC(:,idPsi);
NADH2=C_SUC(:,iNADHm);
UQH2m2=C_SUC(:,iUQH2m);
CytoCm2=C_SUC(:,iCytCred);
Oxygen2=C_SUC(:,iO2);
ADPe2=C_SUC(:,iADPe);
ADPm2=C_SUC(:,iADPm);
R123e_SUC=C_SUC(:,iR123e);
ReBe_SUC=C_SUC(:,iReBe);
R123m_SUC=C_SUC(:,iR123m);
time=0.5:20/60:790/60+0.5;

%---------------------

R123_PM_data=xlsread('R123all.xlsx',1,'A5:E1190');
R123_SUC_data=xlsread('R123all.xlsx',2,'A5:E1190');
%------------

dsfr123=25;    %downsample frequency
T0_R123=1;    %starting time in sec
Nor_PM=mean(mean(R123_PM_data(1:20,:))); % firt 20 points average
Nor_SUC=mean(mean(R123_SUC_data(1:20,:))); % firt 20 points average

R123_PM_data=R123_PM_data/Nor_PM*200e-9;
% R123_PM_data=R123_PM_data*5.404054366078463e-13/2;
R123_PM_mean=mean(R123_PM_data');
R123_PM_se=std(R123_PM_data')/sqrt(5);

R123_SUC_data=R123_SUC_data/Nor_SUC*200e-9;
% R123_SUC_data=R123_SUC_data*5.404054366078463e-13/2;
R123_SUC_mean=mean(R123_SUC_data');
R123_SUC_se=std(R123_SUC_data')/sqrt(5);
%-------------------------------------------------
%Downsample
R123_PM_mean=R123_PM_mean(1:dsfr123:end);
R123_SUC_mean=R123_SUC_mean(1:dsfr123:end);
R123_SUC_se=(R123_SUC_se(1:dsfr123:end));
R123_PM_se=(R123_PM_se(1:dsfr123:end));
R123_time_PM=1:dsfr123:1*dsfr123*length(R123_PM_mean);
R123_time_SUC=1:1*dsfr123:1*dsfr123*length(R123_SUC_mean);



set(figure(1),'Units','inches','Position',[2 2 5 4]) 
set(figure(2),'Units','inches','Position',[1 1 5 4]) 
set(figure(3),'Units','inches','Position',[1 1 5 4]) 
figure(1)
            plot(T+2,Dfi1,'b',T+2,Dfi2,'r','linewidth',2) 
             xlabel('Time (min)')
             ylabel('Membrane Potential (mV)')
      ylim([100 180])     
      xlim([2 20])
legend('PYR+MAL','SUC')
legend boxoff
 title('B: Membrane Potential (Simulation)')
   box off
set(gcf,'color','w')
set(gca,'Fontsize',text_size)


PM_color=[0.3 0.3 1];
SUC_color=[1 0.3 0.3];
              figure(2)
              subplot(2,1,1)
   errorbar(R123_time_PM/60,R123_PM_mean*1e9,R123_PM_se*1e9,'d','marker','*','markersize',mk_size,...
                 'markeredgecolor',PM_color,...
                'color',PM_color,'linewidth',1) 
         ylim([90 220])
         xlim([0 20])
         hold on
           plot(T+2,R123e_PM*1e9,'linewidth',2,'color','k')
          title('A: R123 Concentration')
          box off
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
xlim([0 20])
       hold off
       
           subplot(2,1,2)
           errorbar(R123_time_SUC/60,R123_SUC_mean*1e9,R123_SUC_se*1e9,'d','marker','*','markersize',mk_size,...
                 'markeredgecolor','r',...
                'color','r','linewidth',1) 
    hold on
       plot(T+2,R123e_SUC*1e9,'linewidth',2,'color','k')
               xlabel('Time (min)')
             box off
      set(gcf,'color','w')
set(gca,'Fontsize',text_size)
  
    hold off
             ylim([90 220])
             xlim([0 20])
%-----------------------------------------------------
              figure(3)   
            errorbar(R123_time_PM/60,R123_PM_mean,R123_PM_se,'d','marker','*','markersize',mk_size,...
                 'markeredgecolor',PM_color,...
                'color',PM_color,'linewidth',1) 
           
         hold on
  errorbar(R123_time_SUC/60,R123_SUC_mean,R123_SUC_se,'d','marker','*','markersize',mk_size,...
                 'markeredgecolor',SUC_color,...
                'color',SUC_color,'linewidth',1) 
         hold on
plot(T+2,R123e_PM,'k',T+2,R123e_SUC,'k','linewidth',2)
         xlabel('Time (min)')
         ylabel('Normalized Intensity')
         ylim([0.0 220e-9])
    hold off

%       