
clear all
close all
clc
format long
%%  Parameter Setup
size_index=dlmread('text_size.txt');
text_size=size_index(1);
text_size2=size_index(2);
line_width=size_index(3);
marker_size=size_index(4);
global  Tem C F_con R_con Ve Vm Vi ROTi AAi closed_system...
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
 Tem=298.15-2; %K      23 degree
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 Ve  = 1e-3; %L    
Vm  =1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
ROTi=1;
AAi=1;

closed_system=1;
% Temperature correction
%--------------------------
bufferpH=7.2;
Pi_buffer=5e-3;

%----------------------------
% Temperature correction---------------
% Q10=2.5*ones(25,1);
% Tem_standard=303.15;  %K  30
% Q10_factor=Q10.^((Tem-Tem_standard)/10);
Q10=2.5;
Tem_standard=303.15;  %K  30
Q10_factor=Q10.^((Tem-Tem_standard)/10);

%------------------------------------------

% Q10_factor
Para=ones(27,1);

 Para(1:25)=1*Q10_factor.*Para(1:25);
 
Para(25)=1*0.4*Para(25); %leak is 40% of Evan's 
 UP_factor=7.5;
lengthofS3_factor1=0.90;%10% percent of ADP is consumed, used to calculate peak oxygen consumption rate
lengthofS3_factor2=0.05;%95% percent of ADP is consumed
  %%  Define t_step and t_final
t_step      =   0.01;   %min

%% Run Model
time1=3.7; %Length of State 2  (time before adding ADP)
time2=3;
time3=2;
tic
ADP_add=100e-6;
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:42]);
%options = odeset('NonNegative',[1:42]);
substrates=6;
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iMALe)=5e-3;
IC(iPYRe)=10e-3;
IC(iPie)=Pi_buffer;
% IC(iO2)=293e-6;
Para_PM=Para;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para_PM,2);

IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+ADP_add; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para_PM,2);
IC3=C2(end,:);
 Para_PM(25)= 1*UP_factor*Para_PM(25); %leak
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para_PM,2);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C=[C1; C2(2:end,:);C3(2:end,:)];
Oc_PM=1e6*C(:,iO2);
toc
%-----------------CALCULATE length of S3
index_O2end1=find(C2(:,iADPe)<lengthofS3_factor1*ADP_add,1,'first');
index_O2end2=find(C2(:,iADPe)<lengthofS3_factor2*ADP_add,1,'first');
%------------------------
 O_s2_PM=1e9*(Ve+Vm+Vi)*(C1(end-50,iO2)-C1(end,iO2))/(T1(end)-T1(end-50));%%nmoles/min/mg
 O_s3_PM=1e9*(Ve+Vm+Vi)*(C2(1,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(1));
 O_s4_PM=1e9*(Ve+Vm+Vi)*(C2(index_O2end2+2,iO2)-C2(index_O2end2+7,iO2))/(T2(index_O2end2+7)-T2(index_O2end2+2));


clear C C1 C2

substrates=3;
Para_SUC=Para;
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iSUCe)=7e-3;
IC(iPie)=Pi_buffer;
% IC(iO2)=293e-6;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para_SUC,2);

IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)=ADP_add; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para_SUC,2);
IC3=C2(end,:);
 Para_SUC(25)= UP_factor*Para_SUC(25); %leak
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para_SUC,2);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C=[C1; C2(2:end,:);C3(2:end,:)];
Oc_SUC=1e6*C(:,iO2);

%-----------------CALCULATE length of S3
index_O2end1=find(C2(:,iADPe)<lengthofS3_factor1*ADP_add,1,'first');
index_O2end2=find(C2(:,iADPe)<lengthofS3_factor2*ADP_add,1,'first');
 O_s2_SUC=1e9*(Ve+Vm+Vi)*(C1(end-50,iO2)-C1(end,iO2))/(T1(end)-T1(end-50));%%nmoles/min/mg
 O_s3_SUC=1e9*(Ve+Vm+Vi)*(C2(1,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(1));
 O_s4_SUC=1e9*(Ve+Vm+Vi)*(C2(index_O2end2+20,iO2)-C2(index_O2end2+30,iO2))/(T2(index_O2end2+30)-T2(index_O2end2+20));

  %----------------------------------------------------
substrates=5;
Para_GM=Para;
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iPie)=Pi_buffer;

[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para_GM,2);

IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= ADP_add; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para_GM,2);
IC3=C2(end,:); 
 Para_GM(25)= UP_factor*Para_GM(25); %leak
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para_GM,3);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C=[C1; C2(2:end,:);C3(2:end,:)];
Oc_GM=1e6*C(:,iO2);
index_O2end1=find(C2(:,iADPe)<lengthofS3_factor1*ADP_add,1,'first');
index_O2end2=find(C2(:,iADPe)<lengthofS3_factor2*ADP_add,1,'first');
 O_s2_GM=1e9*(Ve+Vm+Vi)*(C1(end-50,iO2)-C1(end,iO2))/(T1(end)-T1(end-50));%%nmoles/min/mg
 O_s3_GM=1e9*(Ve+Vm+Vi)*(C2(1,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(1));
 O_s4_GM=1e9*(Ve+Vm+Vi)*(C2(index_O2end2+10,iO2)-C2(index_O2end2+20,iO2))/(T2(index_O2end2+20)-T2(index_O2end2+10));

set(figure(5),'Units','inches','Position',[12 2 5 4]) 
set(figure(6),'Units','inches','Position',[12 2 5 4]) 
set(figure(9),'Units','inches','Position',[12 2 5 4]) 

set(figure(11),'Units','inches','Position',[12 2 5 4]) 

%  %-----------------------------------------------
oxi_N8=xlsread('O2_results_n8.xlsx',1);
fs=8;   %downsample rate
MEAN_PM=oxi_N8(2:fs:end,1);
SE_PM=oxi_N8(2:fs:end,2);% 
MEAN_SUC=oxi_N8(2:fs:end,3);
SE_SUC=oxi_N8(2:fs:end,4);
MEAN_GM=oxi_N8(8:fs:end,6);
  Ocsmp_PM=[O_s2_PM,O_s3_PM,O_s4_PM];
  Ocsmp_SUC=[O_s2_SUC,O_s3_SUC,O_s4_SUC];
  Ocsmp_GM=[O_s2_GM,O_s3_GM,O_s4_GM];
    oxi_SUC=[5.305 30.17 5.68];
 SUC_error=[0.9 3.48 0.9];
   oxi_PM=[3.19 17.73 3.62];
 PM_error=[0.08 2.25 0.08];
    oxi_GM=[2.25 21.15 3.8];
 GM_error=[0 0 0];

   %------------------------------------
     figure(5)  %PM figure
   GM_group =  [ oxi_GM' Ocsmp_GM']
   GM_group_SE=[GM_error' [0 0 0]';]
   GM_index=[0.85 1.85 2.85];
   h1=bar(GM_group)
   box off
   ylabel('OCR (nmol/min/mg)')
   title('MAL+GLU')
   set(gcf,'color','w')
   set(gca,'Fontsize',text_size)
set(gca,'XTickLabel',{'State 2','State 3', 'State 4',})
set(h1(1),'FaceColor','b')%--SET Colors
set(h1(2),'FaceColor','r')
set(h1,'EdgeColor','none')
   hold on 
   errorbar(GM_index,oxi_GM,GM_error,'d','marker','none','markersize',5,...
                 'markeredgecolor','b',...
                'color','b','linewidth',2,'CapSize',10)
              legend(h1,'Data','Model')
              ylim([0 40])
               legend boxoff
            hold off

            
                 
time_N8=1:3*fs:(3*fs*length(MEAN_PM));
time_N8=time_N8/60;

      %----------group bar plot  
      set(figure(7),'Units','inches','Position',[12 2 5 4]) 
set(figure(8),'Units','inches','Position',[12 2 5 4]) 
      figure(7)
   PM_group =  [ oxi_PM' Ocsmp_PM']
   PM_group_SE=[PM_error' [0 0 0]';]
   PM_index=[0.85 1.85 2.85];
   h1=bar(PM_group)
   box off
   ylabel('OCR (nmol/min/mg)')
   title('PYR+MAL')
   set(gcf,'color','w')
   set(gca,'Fontsize',text_size)
set(gca,'XTickLabel',{'State 2','State 3', 'State 4',})
set(h1(1),'FaceColor','b')%--SET Colors
set(h1(2),'FaceColor','r')
set(h1,'EdgeColor','none')
   hold on 
   errorbar(PM_index,oxi_PM,PM_error,'d','marker','none','markersize',5,...
                 'markeredgecolor','b',...
                'color','b','linewidth',2,'CapSize',10)
              legend(h1,'Data','Model')
              ylim([0 40])
               legend boxoff
            hold off
            
             figure(11)
   SUC_group =  [ oxi_SUC' Ocsmp_SUC']
   SUC_group_SE=[SUC_error' [0 0 0]';]
   SUC_index=[0.85 1.85 2.85];
   h1=bar(SUC_group)
   
   box off
   ylabel('OCR (nmol/min/mg)')
   title('SUC')
   set(gcf,'color','w')
   set(gca,'Fontsize',text_size)
set(gca,'XTickLabel',{'State 2','State 3', 'State 4',})
set(h1(1),'FaceColor','b')%--SET Colors
set(h1(2),'FaceColor','r')
set(h1,'EdgeColor','none')
   hold on 
   errorbar(SUC_index,oxi_SUC,SUC_error,'d','marker','none','markersize',5,...
                 'markeredgecolor','b',...
                'color','b','linewidth',2,'CapSize',10)
             ylim([0 40])
              legend(h1,'Data','Model')
               legend boxoff
            hold off
            
  time_GM=1:3*fs:(3*fs*length(MEAN_GM));
time_GM=time_GM/60;


figure(8) 

h1=plot(T,Oc_PM,'r','LineWidth',2)
box off
title('PYR+MAL')

xlabel('Time (min)')
ylabel('Oxygen Concentration (\muM)')
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
xlim([0 9])
ylim([180 300])
hold on
h2=errorbar(time_N8-2,MEAN_PM,SE_PM,'d','marker','square','markersize',10,...
                 'markeredgecolor','b',...
                'color','b','linewidth',2)

   legend([h2 h1],'Data','Model')
   legend boxoff
        
figure(9) 

h1=plot(T,Oc_SUC,'r','LineWidth',2)
box off
title('SUC')

xlabel('Time (min)')
ylabel('Oxygen Concentration (\muM)')
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
xlim([0 9])
ylim([180 300])
hold on
h2= errorbar(time_N8-2,MEAN_SUC,SE_SUC,'d','marker','square','markersize',10,...
                 'markeredgecolor','b',...
                'color','b','linewidth',2)
       legend([h2 h1],'Data','Model')
       legend boxoff
        hold off
                  
  set(figure(10),'Units','inches','Position',[12 2 5 4])           
      
figure(10) 

plot(time_GM,MEAN_GM,'bs',T,Oc_GM,'r','LineWidth',line_width,'markersize',10)
box off
title('MAL+GLU')
legend('Data','Model')

xlabel('Time (min)')
ylabel('Oxygen Concentration (\muM)')
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
xlim([0 9])
ylim([180 300])
            legend boxoff
            
    