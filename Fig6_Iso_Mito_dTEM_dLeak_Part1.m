
clear all
clc
format long

%----
text_size=18;
mk_size=12;
line_width=1.5;
%----------------------
%%  Parameter Setup
bufferpH=7.2;
global  Tem C F_con R_con Ve Vm Vi ROTi AAi closed_system...
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

R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 Ve  =   1e-3; %L    
Vm  = 1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
ROTi=1;
AAi=1;
closed_system=1;% closed system

  Tem_standard=303.15;  %K  30
substrates=6;
lengthofS3_factor1=0.90;%10% percent of ADP is consumed
lengthofS3_factor2=0.80;%95% percent of ADP is consumed
  %%  Define t_step and t_final
t_step      =   0.02;   %min
%% Run Simulation
 dTEM=[(303.15-7):2:(303.15+7)];
 dLeak=[0.4:0.3:1.6];

 O_s3=zeros(5,1);
 O_s2=zeros(5,1);
 Para=ones(27,1);
for i=1:1:length(dTEM)
    for j=1:1:length(dLeak)
 ADP_amount=300e-6;

 Tem=dTEM(i);
 leak_factor=dLeak(j);
Q10=2.5*ones(25,1);
Tem_standard=303.15;  %K  30
Q10_factor=Q10.^((Tem-Tem_standard)/10);
%------------------------------------------
 dPara=Para;
dPara(1:25)=1*Q10_factor.*dPara(1:25);
dPara(25)=leak_factor*dPara(25);
% Tem=310.15;
time1=3; %Length of State 2  (time before adding ADP)
time2=7;
time3=2;

IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iPYRe)=10e-3;
IC(iMALe)=5e-3;

options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:44]);

[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,dPara,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)=IC2(iADPe)+ ADP_amount; %add ADP, Unit(Molar) 

[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,dPara,2);


T=[T1; T2(2:end)+time1;];

C=[C1; C2(2:end,:);];

index_O2end1=find(C2(:,iADPe)<lengthofS3_factor1*ADP_amount,1,'first') %10% oxygen consumed

 O_s2(i,j)=1e9*(Ve+Vm+Vi)*(C1(time1/t_step/2,iO2)-C1(time1/t_step,iO2))/(T1(time1/t_step)-T1(time1/t_step/2));%nmols
 O_s3(i,j)=1e9*(Ve+Vm+Vi)*(C2(2,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(2));

clear C dPara
    end
end

%%----------------------
set(figure(1),'Units','inches','Position',[0.2 0.1 5 4]) 
set(figure(2),'Units','inches','Position',[5.2 0.1 5 4]) 
set(figure(3),'Units','inches','Position',[0.2 4.1 5 4]) 
set(figure(4),'Units','inches','Position',[5.2 4.1 5 4]) 
set(figure(5),'Units','inches','Position',[10.2 4.1 5 4]) 
set(figure(6),'Units','inches','Position',[10.2 4.1 5 4]) 
%------------------------------------
 meanOCR_SUC_23=[5.305 30.17]; %nmoles/min/mg
 seOCR_SUC_23=[0.9 3.48];
 meanOCR_SUC_30=[26.3 68.9]; %nmoles/min/mg
 seOCR_SUC_30=[4 4];
  meanOCR_PYRMAL_23=[3.19 17.73];
 seOCR_PYRMAL_23=[0.08 2.25];
 meanOCR_PYRMAL_30=[12.6 35];
 seOCR_PYRMAL_30=[1.7 2.4];
 mean_SUCROT=[56.11 77.63 56.11];
 se_SUCROT=[6.7 10.9 6.7];
 
figure(1)
h1=plot(dTEM-273.15,O_s2(:,1),'-.r',dTEM-273.15,O_s2(:,2),'-.g',dTEM-273.15,O_s2(:,3),'-.b',dTEM-273.15,O_s2(:,4),'-.m',dTEM-273.15,O_s2(:,5),'-.c','LineWidth',line_width,'MarkerSize',8)
hold on
h2=plot(dTEM-273.15,O_s3(:,1),'r',dTEM-273.15,O_s3(:,2),'g',dTEM-273.15,O_s3(:,3),'b',dTEM-273.15,O_s3(:,4),'m',dTEM-273.15,O_s3(:,5),'c','LineWidth',line_width,'MarkerSize',8)
box off
xlim([20 40])
ylim([0 130])
hold off
xlabel(['Temperatue ' char(176) 'C'])
ylabel('OCR (nmol/min/mg)')
set(gcf,'color','w')
set(gca,'Fontsize',text_size)

if substrates==3
title('States 2 and 3 OCR (SUC)','FontSize',text_size)
hold on
 h4=errorbar([30 30],meanOCR_SUC_30,seOCR_SUC_30,'d','marker','square','markersize',8,...
                 'markeredgecolor','k',...
                'color','k','linewidth',2)
            hold on
 h4=errorbar([23 23],meanOCR_SUC_23,seOCR_SUC_23,'d','marker','square','markersize',8,...
                 'markeredgecolor','k',...
                'color','k','linewidth',2)
            hold off
            legend(h2,'40% Leak','70% Leak','100% Leak','130% Leak','170% Leak')
elseif substrates==6
    title('States 2 and 3 OCR (PYR+MAL)','FontSize',text_size)
hold on
 h4=errorbar([30 30],meanOCR_PYRMAL_30,seOCR_PYRMAL_30,'d','marker','square','markersize',8,...
                 'markeredgecolor','k',...
                'color','k','linewidth',2)
            hold on
 h4=errorbar([23 23],meanOCR_PYRMAL_23,seOCR_PYRMAL_23,'d','marker','square','markersize',8,...
                 'markeredgecolor','k',...
                'color','k','linewidth',2)
            hold off
            legend(h2,'40% Leak','70% Leak','100% Leak','130% Leak','170% Leak')
end
legend boxoff
figure(2)
plot(dTEM-273.15,O_s3(:,1)./O_s2(:,1),'r',dTEM-273.15,O_s3(:,2)./O_s2(:,2),'g',dTEM-273.15,O_s3(:,3)./O_s2(:,3),'b',dTEM-273.15,O_s3(:,4)./O_s2(:,4),'m',dTEM-273.15,O_s3(:,5)./O_s2(:,5),'c','LineWidth',line_width,'MarkerSize',8)
box off
xlabel(['Temperatue ' char(176) 'C'])
ylabel('RCI')
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
legend('40% Leak','70% Leak','100% Leak','130% Leak','170% Leak')
legend boxoff

if substrates==3
title('RCI (SUC)','FontSize',text_size)
elseif substrates==6
    title('RCI (PYR+MAL)','FontSize',text_size)
end
