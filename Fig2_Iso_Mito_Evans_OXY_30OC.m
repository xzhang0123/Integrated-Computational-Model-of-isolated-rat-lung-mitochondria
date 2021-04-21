
clear all
close all
clc
format long
%--------------
size_index=dlmread('text_size.txt');
text_size=size_index(1);
text_size2=size_index(2);
line_width=size_index(3);
marker_size=size_index(4);
shade_area=[0.75 0.75 0.75];
%--------------------------------
 
 %--------EXPERIMENT DATA-----------
 meanOCR_SUC=[26.3 68.9]; %nmoles/min/mg
 seOCR_SUC=[4 4];
 meanOCR_PYRMAL=[12.6 35];
 seOCR_PYRMAL=[1.7 2.4];
 meanOCR_GM=[22.8 52.15 41.62/2];
 seOCR_GM=[3 6.4 5.33/2];
 mean_SUCROT=[56.11 77.63 56.11];
 se_SUCROT=[6.7 10.9 6.7];
 meanOCR_AKGPYR=[10.4 30.1];
 seOCR_AKGPYR=[1.4 5.3];
  meanOCR_CITPYR=[12.8 42.3];
 seOCR_CITPYR=[2.2 7.2];
 %-------------Model---------------------------
%%  Parameter Setup
substrates=1;
bufferpH=7.4;
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
 Tem=303.15; %K      30 oC
 %Tem=298.15; %K      25 oC   25 degree
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
Ve  =   1.9/0.4*1e-3; %L   
Vm  = 1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
ROTi=1;%CI inhibited by ROT  percentage
AAi=1;
closed_system=1;
lengthofS3_factor1=0.99;%1% percent of ADP is consumed
lengthofS3_factor2=0.90;%10% percent of ADP is consumed
index_O2start1=2;
% Temperature correction---------------
Q10=2.5*ones(25,1);
Tem_standard=303.15;  %K  30
Q10_factor=Q10.^((Tem-Tem_standard)/10);
%------------------------------------------
Para=ones(27,1);


 Para(1:25)=1*Q10_factor.*Para(1:25);




  %%  Model settings
t_step      =   0.01;   %min
time1=1; %Length of State 2  (time before adding ADP)
time2=0.5;
time3=0.5;
tic
ADP_amount=300e-6;
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:42]);
%% CITPYR Model
IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iCITe)=10e-3;
IC(iPYRe)=5e-3;
IC(iPie)=4e-3;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+1*ADP_amount; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
IC3=C2(end,:);
IC3(iADPe)= IC3(iADPe)+0*ADP_amount; %add ADP, Unit(Molar) 
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para,2);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C=[C1; C2(2:end,:);C3(2:end,:)];
OC_CITPYR=C(:,iO2);
%-----------------CALCULATE length of S3
index_O2end1=find(C2(:,iADPe)<lengthofS3_factor1*ADP_amount,1,'first');
index_O2end2=find(C2(:,iADPe)<lengthofS3_factor2*ADP_amount,1,'first');
 O_s2_CITPYR=1e9*(Ve+Vm+Vi)*(C1(end-50,iO2)-C1(end,iO2))/(T1(end)-T1(end-50));%nmols
 O_s3_CITPYR=1e9*(Ve+Vm+Vi)*(C2(index_O2start1,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(index_O2start1));
  OCR_SIM_CITPYR=[O_s2_CITPYR,O_s3_CITPYR];
 %%---------------------------SUC Model
 %% 
 IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iSUCe)=5e-3;
IC(iPYRe)=5e-3; %PYR also added into buffer, but no significant diference, Evans et al. 
IC(iPie)=4e-3;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+1*ADP_amount; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
IC3=C2(end,:);
IC3(iADPe)= IC3(iADPe)+0*ADP_amount; %add ADP, Unit(Molar) 
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para,2);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C=[C1; C2(2:end,:);C3(2:end,:)];
OC_SUC=C(:,iO2);
%-----------------CALCULATE length of S3
index_O2end1=find(C2(:,iADPe)<lengthofS3_factor1*ADP_amount,1,'first');
index_O2end2=find(C2(:,iADPe)<lengthofS3_factor2*ADP_amount,1,'first');
 O_s2_SUC=1e9*(Ve+Vm+Vi)*(C1(end-50,iO2)-C1(end,iO2))/(T1(end)-T1(end-50));%nmols
 O_s3_SUC=1e9*(Ve+Vm+Vi)*(C2(index_O2start1,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(index_O2start1));
  OCR_SIM_SUC=[O_s2_SUC,O_s3_SUC];
  %--------------------------------------
  %%  %%---------------------------aKG PYR Model
 %% 
 IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iPYRe)=5e-3;
IC(iaKGe)=5e-3;
IC(iPie)=4e-3;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+1*ADP_amount; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
IC3=C2(end,:);
IC3(iADPe)= IC3(iADPe)+0*ADP_amount; %add ADP, Unit(Molar) 
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para,2);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C=[C1; C2(2:end,:);C3(2:end,:)];
OC_AKGPYR=C(:,iO2);
%-----------------CALCULATE length of S3
index_O2end1=find(C2(:,iADPe)<lengthofS3_factor1*ADP_amount,1,'first');
index_O2end2=find(C2(:,iADPe)<lengthofS3_factor2*ADP_amount,1,'first');
 O_s2_AKGPYR=1e9*(Ve+Vm+Vi)*(C1(end-50,iO2)-C1(end,iO2))/(T1(end)-T1(end-50));%nmols
 O_s3_AKGPYR=1e9*(Ve+Vm+Vi)*(C2(index_O2start1,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(index_O2start1));
  OCR_SIM_AKGPYR=[O_s2_AKGPYR,O_s3_AKGPYR];
  %--------------------------------------
    %--------------------------------------
  %%  %%---------------------------PYRMAL Model
 %% 
 IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iPYRe)=5e-3;
IC(iMALe)=5e-3;
IC(iPie)=4e-3;
[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,Para,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)= IC2(iADPe)+1*ADP_amount; %add ADP, Unit(Molar) 
[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,Para,2);
IC3=C2(end,:);
IC3(iADPe)= IC3(iADPe)+0*ADP_amount; %add ADP, Unit(Molar) 
[T3,C3]= ode15s(@odeq,[0:t_step:time3],IC3,options,substrates,Para,2);
T=[T1; T2(2:end)+time1;T3(2:end)+time1+time2];
C=[C1; C2(2:end,:);C3(2:end,:)];
OC_PYRMAL=C(:,iO2);
%-----------------CALCULATE length of S3
index_O2end1=find(C2(:,iADPe)<lengthofS3_factor1*ADP_amount,1,'first');
index_O2end2=find(C2(:,iADPe)<lengthofS3_factor2*ADP_amount,1,'first');
 O_s2_PYRMAL=1e9*(Ve+Vm+Vi)*(C1(end-50,iO2)-C1(end,iO2))/(T1(end)-T1(end-50));%nmols
 O_s3_PYRMAL=1e9*(Ve+Vm+Vi)*(C2(index_O2start1,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(index_O2start1));
  OCR_SIM_PYRMAL=[O_s2_PYRMAL,O_s3_PYRMAL];
  %--------------------------------------
set(figure(3),'Units','inches','Position',[2 2 5 4]) 
set(figure(4),'Units','inches','Position',[2 2 5 4]) 
set(figure(5),'Units','inches','Position',[2 2 5 4]) 
set(figure(6),'Units','inches','Position',[2 2 5 4]) 
figure(3)

fill_x1=linspace(0,time1,10)';
fill_y1_f=-(meanOCR_CITPYR(1)-seOCR_CITPYR(1))/(Ve+Vi+Vm)./1e3.*fill_x1+296;  %forward state 2  up boundary
fill_y1_b=-(meanOCR_CITPYR(1)+seOCR_CITPYR(1))/(Ve+Vi+Vm)./1e3.*fill_x1+296;  %reverse state 2
%-----
fill_x2=linspace(time1,time1+1,10)';
fill_y2_f=-(meanOCR_CITPYR(2)-seOCR_CITPYR(2))/(Ve+Vi+Vm)./1e3.*(fill_x2-time1)+...
    -(meanOCR_CITPYR(1)-seOCR_CITPYR(1))/(Ve+Vi+Vm)./1e3.*time1+296; %forward state 3, intercetpt at time1
fill_y2_b=-(meanOCR_CITPYR(2)+seOCR_CITPYR(2))/(Ve+Vi+Vm)./1e3.*(fill_x2-time1)+...
    -(meanOCR_CITPYR(1)+seOCR_CITPYR(1))/(Ve+Vi+Vm)./1e3.*time1+296; %forward state 3, intercetpt at time1
fill([fill_x1;flipud(fill_x1)],[fill_y1_f;flipud(fill_y1_b)],shade_area,'linestyle','none');
 set(gcf,'color','w')
set(gca,'Fontsize',text_size)
legend boxoff 
box off
title('CIT+PYR')
hold on
h2_f3=fill([fill_x2;flipud(fill_x2)],[fill_y2_f;flipud(fill_y2_b)],shade_area,'linestyle','none');
hold on
h1_f3=plot(T,OC_CITPYR*1e6,'k','linewidth',line_width)
xlabel('Time (min)')
ylabel('Oxygen Concentration (\muM)')
ylim([275 296])
 legend([h1_f3,h2_f3],'Model','Data');
 legend boxoff
hold off

%----------------------------------    
figure(4)

fill_x1=linspace(0,time1,10)';
fill_y1_f=-(meanOCR_SUC(1)-seOCR_SUC(1))/(Ve+Vi+Vm)./1e3.*fill_x1+296;  %forward state 2  up boundary (slow)
fill_y1_b=-(meanOCR_SUC(1)+seOCR_SUC(1))/(Ve+Vi+Vm)./1e3.*fill_x1+296;  %reverse state 2
%-----
fill_x2=linspace(time1,time1+1,10)';
fill_y2_f=-(meanOCR_SUC(2)-seOCR_SUC(2))/(Ve+Vi+Vm)./1e3.*(fill_x2-time1)+...
    -(meanOCR_SUC(1)-seOCR_SUC(1))/(Ve+Vi+Vm)./1e3.*time1+296; %forward state 3, intercetpt at time1
fill_y2_b=-(meanOCR_SUC(2)+seOCR_SUC(2))/(Ve+Vi+Vm)./1e3.*(fill_x2-time1)+...
    -(meanOCR_SUC(1)+seOCR_SUC(1))/(Ve+Vi+Vm)./1e3.*time1+296; %forward state 3, intercetpt at time1
fill([fill_x1;flipud(fill_x1)],[fill_y1_f;flipud(fill_y1_b)],shade_area,'linestyle','none');
 set(gcf,'color','w')
set(gca,'Fontsize',text_size)
legend boxoff 
box off
title('SUC')
hold on
h2_f4=fill([fill_x2;flipud(fill_x2)],[fill_y2_f;flipud(fill_y2_b)],shade_area,'linestyle','none');
hold on
h1_f4=plot(T,OC_SUC*1e6,'k','linewidth',line_width)
xlabel('Time (min)')
ylabel('Oxygen Concentration (\muM)')
ylim([275 296])
 legend([h1_f4,h2_f4],'Model','Data');
 legend boxoff
hold off
%----------------------------------    
figure(5)


fill_x1=linspace(0,time1,10)';
fill_y1_f=-(meanOCR_PYRMAL(1)-seOCR_PYRMAL(1))/(Ve+Vi+Vm)./1e3.*fill_x1+296;  %forward state 2  up boundary (slow)
fill_y1_b=-(meanOCR_PYRMAL(1)+seOCR_PYRMAL(1))/(Ve+Vi+Vm)./1e3.*fill_x1+296;  %reverse state 2
%-----
fill_x2=linspace(time1,time1+1,10)';
fill_y2_f=-(meanOCR_PYRMAL(2)-seOCR_PYRMAL(2))/(Ve+Vi+Vm)./1e3.*(fill_x2-time1)+...
    -(meanOCR_PYRMAL(1)-seOCR_PYRMAL(1))/(Ve+Vi+Vm)./1e3.*time1+296; %forward state 3, intercetpt at time1
fill_y2_b=-(meanOCR_PYRMAL(2)+seOCR_PYRMAL(2))/(Ve+Vi+Vm)./1e3.*(fill_x2-time1)+...
    -(meanOCR_PYRMAL(1)+seOCR_PYRMAL(1))/(Ve+Vi+Vm)./1e3.*time1+296; %forward state 3, intercetpt at time1
fill([fill_x1;flipud(fill_x1)],[fill_y1_f;flipud(fill_y1_b)],shade_area,'linestyle','none');
title('PYR+MAL')
 set(gcf,'color','w')
set(gca,'Fontsize',text_size)
legend boxoff 
box off
hold on
h2_f5=fill([fill_x2;flipud(fill_x2)],[fill_y2_f;flipud(fill_y2_b)],shade_area,'linestyle','none');
hold on
h1_f5=plot(T,OC_PYRMAL*1e6,'k','linewidth',line_width)
xlabel('Time (min)')
ylabel('Oxygen Concentration (\muM)')
ylim([275 296])
 legend([h1_f5,h2_f5],'Model','Data');
 legend boxoff
hold off
%----------------------------------    
figure(6)

fill_x1=linspace(0,time1,10)';
fill_y1_f=-(meanOCR_AKGPYR(1)-seOCR_AKGPYR(1))/(Ve+Vi+Vm)./1e3.*fill_x1+296;  %forward state 2  up boundary (slow)
fill_y1_b=-(meanOCR_AKGPYR(1)+seOCR_AKGPYR(1))/(Ve+Vi+Vm)./1e3.*fill_x1+296;  %reverse state 2
%-----
fill_x2=linspace(time1,time1+1,10)';
fill_y2_f=-(meanOCR_AKGPYR(2)-seOCR_AKGPYR(2))/(Ve+Vi+Vm)./1e3.*(fill_x2-time1)+...
    -(meanOCR_AKGPYR(1)-seOCR_AKGPYR(1))/(Ve+Vi+Vm)./1e3.*time1+296; %forward state 3, intercetpt at time1
fill_y2_b=-(meanOCR_AKGPYR(2)+seOCR_AKGPYR(2))/(Ve+Vi+Vm)./1e3.*(fill_x2-time1)+...
    -(meanOCR_AKGPYR(1)+seOCR_AKGPYR(1))/(Ve+Vi+Vm)./1e3.*time1+296; %forward state 3, intercetpt at time1
fill([fill_x1;flipud(fill_x1)],[fill_y1_f;flipud(fill_y1_b)],shade_area,'linestyle','none');
title('PYR+AKG')
 set(gcf,'color','w')
set(gca,'Fontsize',text_size)
hold on
h2_f6=fill([fill_x2;flipud(fill_x2)],[fill_y2_f;flipud(fill_y2_b)],shade_area,'linestyle','none');
hold on
h1_f6=plot(T,OC_AKGPYR*1e6,'k','linewidth',line_width)
xlabel('Time (min)')
ylabel('Oxygen Concentration (\muM)')
legend boxoff 
box off
ylim([275 296])
hold off
 legend([h1_f6,h2_f6],'Model','Data');
 legend boxoff
%---------------------------------

set(figure(11),'Units','inches','Position',[2 2 6 4.5]) 
set(figure(12),'Units','inches','Position',[2 2 6 4.5]) 
set(figure(13),'Units','inches','Position',[2 2 6 4.5]) 
set(figure(14),'Units','inches','Position',[2 2 6 4.5]) 
 figure(11)
   CITPYR_group =  [ meanOCR_CITPYR' OCR_SIM_CITPYR']
   CITPYR_group_SE=[seOCR_CITPYR' [0 0]';]
   CITPYR_index=[0.85 1.85];
   h1=bar(CITPYR_group)
   ylim([0 80])
   box off
   ylabel('OCR (nmol/min/mg)')
  % title('CIT+PYR')
   set(gcf,'color','w')
   set(gca,'Fontsize',1.5*text_size)
set(gca,'XTickLabel',{'State 2','State 3'})
set(h1(1),'FaceColor','b')%--SET Colors
set(h1(2),'FaceColor','r')
set(h1,'EdgeColor','none')
   hold on 
   errorbar(CITPYR_index,meanOCR_CITPYR,seOCR_CITPYR,'d','marker','none','markersize',5,...
                 'markeredgecolor','b',...
                'color','b','linewidth',5,'CapSize',20)
              legend(h1,'Data','Model')
               legend boxoff
            hold off


             figure(12)
   SUC_group =  [ meanOCR_SUC' OCR_SIM_SUC']
   SUC_group_SE=[seOCR_SUC' [0 0]';]
   SUC_index=[0.85 1.85];
   h1=bar(SUC_group)
   ylim([0 80])
   box off
   ylabel('OCR (nmol/min/mg)')
  % title('SUC')
   set(gcf,'color','w')
   set(gca,'Fontsize',1.5*text_size)
set(gca,'XTickLabel',{'State 2','State 3'})
set(h1(1),'FaceColor','b')%--SET Colors
set(h1(2),'FaceColor','r')
set(h1,'EdgeColor','none')
   hold on 
   errorbar(SUC_index,meanOCR_SUC,seOCR_SUC,'d','marker','none','markersize',5,...
                 'markeredgecolor','b',...
                'color','b','linewidth',5,'CapSize',20)
              legend(h1,'Data','Model')
               legend boxoff
            hold off


             figure(13)
   PYRMAL_group =  [ meanOCR_PYRMAL' OCR_SIM_PYRMAL'];
   PYRMAL_group_SE=[seOCR_PYRMAL' [0 0]';];
   PYRMAL_index=[0.85 1.85];
   h1=bar(PYRMAL_group)
   ylim([0 80])
   box off
   ylabel('OCR (nmol/min/mg)')
  % title('PYRMAL')
   set(gcf,'color','w')
   set(gca,'Fontsize',1.5*text_size)
set(gca,'XTickLabel',{'State 2','State 3'})
set(h1(1),'FaceColor','b')%--SET Colors
set(h1(2),'FaceColor','r')
set(h1,'EdgeColor','none')
   hold on 
   errorbar(PYRMAL_index,meanOCR_PYRMAL,seOCR_PYRMAL,'d','marker','none','markersize',5,...
                 'markeredgecolor','b',...
                'color','b','linewidth',5,'CapSize',20)
              legend(h1,'Data','Model')
               legend boxoff
            hold off


   figure(14)
   AKGPYR_group =  [ meanOCR_AKGPYR' OCR_SIM_AKGPYR']
   AKGPYR_group_SE=[seOCR_AKGPYR' [0 0]';]
   AKGPYR_index=[0.85 1.85];
   h1=bar(AKGPYR_group)
   ylim([0 80])
   box off
   ylabel('OCR (nmol/min/mg)')
  % title('PYR+AKG')
   set(gcf,'color','w')
   set(gca,'Fontsize',1.5*text_size)
set(gca,'XTickLabel',{'State 2','State 3'})
set(h1(1),'FaceColor','b')%--SET Colors
set(h1(2),'FaceColor','r')
set(h1,'EdgeColor','none')
   hold on 
   errorbar(AKGPYR_index,meanOCR_AKGPYR,seOCR_AKGPYR,'d','marker','none','markersize',5,...
                 'markeredgecolor','b',...
                'color','b','linewidth',5,'CapSize',20)
              legend(h1,'Data','Model')
               legend boxoff
            hold off

