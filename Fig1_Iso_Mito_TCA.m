
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
mk_size=10;
line_width=1.2;
%-------------------------
%%  Parameter Setup
substrates=6; %PYR+MAL
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
Ve  =   1/0.7*1e-3; %L normalized to 1mg mitochondria
Vm  =   1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
%---------------------
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
   Ti0=5; %time before measuring
time1=10; %Length of State 2  (time before adding ADP)
tic
PYRstep=[0 0.5e-3 1e-3 5e-3 10e-3];
MALstep=[0 0.5e-3 1e-3 5e-3 10e-3];

%---------------
CIT_rel1=zeros(time1/t_step+1,5);
PYR_upt1=zeros(time1/t_step+1,5);
MAL_upt1=zeros(time1/t_step+1,5);

%----Simulation with differnet [PYRe]-----
for i=1:1:5
    IC=Set_Initial_Concentrations(substrates,bufferpH);
    IC(iPYRe)=PYRstep(i);
    IC(iMALe)=5e-3;

options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on');
%       options = odeset('NonNegative',[1:44]);
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
CIT_rel1(:,i)=(C(:,iCITe)-C(1,iCITe))*Ve*1e9;
PYR_upt1(:,i)=1e9*(C1(1,iPYRe)-C1(:,iPYRe))*Ve;
MAL_upt1(:,i)=(C(1,iMALe)-C(:,iMALe))*Ve*1e9;
%-----------------------
clear T1 C1 C T00 C00 T
end
%----------Simulation with differnet [MALe]----------
IC=Set_Initial_Concentrations(substrates,bufferpH);
for i=1:1:5
    
   IC(iMALe)=MALstep(i);
   IC(iPYRe)=5e-3;
[T00,C00] = ode15s(@odeq,[0:t_step:Ti0],IC,options,substrates,Para,2);
IC1=C00(end,:);    
      [T1,C1] = ode15s(@odeq,[0:t_step:time1],IC1,options,substrates,Para,2);
C=C1;
T=T1;
CIT_rel2(:,i)=(C(:,iCITe)-C(1,iCITe))*Ve*1e9;
PYR_upt2(:,i)=1e9*(C(1,iPYRe)-C(:,iPYRe))*Ve;
MAL_upt2(:,i)=(C(1,iMALe)-C(:,iMALe))*Ve*1e9;
clear T1 C1 C
end

toc
clear t t_step t_final

Ax1=xlsread('TCA1975.xlsx','A3:A7');
Ay1=xlsread('TCA1975.xlsx','B3:B7');
Ax2=xlsread('TCA1975.xlsx','A9:A13');
Ay2=xlsread('TCA1975.xlsx','B9:B13');
Ax3=xlsread('TCA1975.xlsx','A15:A19');
Ay3=xlsread('TCA1975.xlsx','B15:B19');
Ax4=xlsread('TCA1975.xlsx','A21:A25');
Ay4=xlsread('TCA1975.xlsx','B21:B25');
Ax5=xlsread('TCA1975.xlsx','A27:A31');
Ay5=xlsread('TCA1975.xlsx','B27:B31');

Bx1=xlsread('TCA1975.xlsx','D3:D7');
By1=xlsread('TCA1975.xlsx','E3:E7');
Bx2=xlsread('TCA1975.xlsx','D9:D13');
By2=xlsread('TCA1975.xlsx','E9:E13');
Bx3=xlsread('TCA1975.xlsx','D15:D19');
By3=xlsread('TCA1975.xlsx','E15:E19');
Bx4=xlsread('TCA1975.xlsx','D21:D25');
By4=xlsread('TCA1975.xlsx','E21:E25');
Bx5=xlsread('TCA1975.xlsx','D27:D31');
By5=xlsread('TCA1975.xlsx','E27:E31');

Cx1=xlsread('TCA1975.xlsx','G3:G7');
Cy1=xlsread('TCA1975.xlsx','H3:H7');
Cx2=xlsread('TCA1975.xlsx','G9:G13');
Cy2=xlsread('TCA1975.xlsx','H9:H13');
Cx3=xlsread('TCA1975.xlsx','G15:G19');
Cy3=xlsread('TCA1975.xlsx','H15:H19');
Cx4=xlsread('TCA1975.xlsx','G21:G25');
Cy4=xlsread('TCA1975.xlsx','H21:H25');
Cx5=xlsread('TCA1975.xlsx','G27:G31');
Cy5=xlsread('TCA1975.xlsx','H27:H31');

Dx1=xlsread('TCA1975.xlsx','J3:J7');
Dy1=xlsread('TCA1975.xlsx','K3:K7');
Dx2=xlsread('TCA1975.xlsx','J9:J13');
Dy2=xlsread('TCA1975.xlsx','K9:K13');
Dx3=xlsread('TCA1975.xlsx','J15:J19');
Dy3=xlsread('TCA1975.xlsx','K15:K19');
Dx4=xlsread('TCA1975.xlsx','J21:J25');
Dy4=xlsread('TCA1975.xlsx','K21:K25');
Dx5=xlsread('TCA1975.xlsx','J27:J31');
Dy5=xlsread('TCA1975.xlsx','K27:K31');

Ex1=xlsread('TCA1975.xlsx','L3:L7');
Ey1=xlsread('TCA1975.xlsx','M3:M7');
Ex2=xlsread('TCA1975.xlsx','L9:L13');
Ey2=xlsread('TCA1975.xlsx','M9:M13');
Ex3=xlsread('TCA1975.xlsx','L15:L19');
Ey3=xlsread('TCA1975.xlsx','M15:M19');
Ex4=xlsread('TCA1975.xlsx','L21:L25');
Ey4=xlsread('TCA1975.xlsx','M21:M25');
Ex5=xlsread('TCA1975.xlsx','L27:L31');
Ey5=xlsread('TCA1975.xlsx','M27:M31');

Fx1=xlsread('TCA1975.xlsx','N3:N7');
Fy1=xlsread('TCA1975.xlsx','O3:O7');
Fx2=xlsread('TCA1975.xlsx','N9:N13');
Fy2=xlsread('TCA1975.xlsx','O9:O13');
Fx3=xlsread('TCA1975.xlsx','N15:N19');
Fy3=xlsread('TCA1975.xlsx','O15:O19');
Fx4=xlsread('TCA1975.xlsx','N21:N25');
Fy4=xlsread('TCA1975.xlsx','O21:O25');
Fx5=xlsread('TCA1975.xlsx','N27:N31');
Fy5=xlsread('TCA1975.xlsx','O27:O31');

set(figure(1),'Units','inches','Position',[2 1 15 4]) 
figure(1)
hold on
subplot(1,3,1)
plot(Bx1,By1,'ro',Bx2,By2,'g*',Bx3,By3,'b>',Bx4,By4,'m+',Bx5,By5,'c^',T,PYR_upt1(:,1),'r',T,PYR_upt1(:,2),'g',T,PYR_upt1(:,3),'b',T,PYR_upt1(:,4),'m',T,PYR_upt1(:,5),'c','LineWidth',2,'MarkerSize',mk_size)
  box off 
  set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
set(gca,'XTick',(0:2:10),'YTick',(0:30:150));
text(7,10,'[MAL]=5mM','FontSize',text_size);
legend('[PYR]=0 mM','[PYR]=0.5 mM','[PYR]=1 mM','[PYR]=5 mM','[PYR]=10 mM','Location','NorthWest')
legend boxoff
title('PYR Uptake')
ylabel('Mass (nmoles/mg)')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',text_size)

subplot(1,3,2)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
plot(Cx1,Cy1,'ro',Cx2,Cy2,'g*',Cx3,Cy3,'b>',Cx4,Cy4,'m+',Cx5,Cy5,'c^',T,MAL_upt1(:,1),'r',T,MAL_upt1(:,2),'g',T,MAL_upt1(:,3),'b',T,MAL_upt1(:,4),'m',T,MAL_upt1(:,5),'c','LineWidth',2,'MarkerSize',mk_size)
 box off  
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
set(gca,'XTick',(0:2:10),'YTick',(0:30:150));
ylim([0 150])
title('MAL Uptake')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',text_size)
subplot(1,3,3)


plot(Ax1,Ay1,'ro',Ax2,Ay2,'g*',Ax3,Ay3,'b>',Ax4,Ay4,'m+',Ax5,Ay5,'c^',T,CIT_rel1(:,1),'r',T,CIT_rel1(:,2),'g',T,CIT_rel1(:,3),'b',T,CIT_rel1(:,4),'m',T,CIT_rel1(:,5),'c','LineWidth',2,'MarkerSize',mk_size)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
set(gca,'XTick',(0:2:10),'YTick',(0:30:150));
ylim([0 150])
box off 
title('CIT Production')
xlim([0 11])
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',text_size)
%----------------------------------------
set(figure(2),'Units','inches','Position',[3 3 15 4]) 
figure(2)
box off 
subplot(1,3,1)
plot(Ex1,Ey1,'ro',Ex2,Ey2,'g*',Ex3,Ey3,'b>',Ex4,Ey4,'m+',Ex5,Ey5,'c^',T,PYR_upt2(:,1),'r',T,PYR_upt2(:,2),'g',T,PYR_upt2(:,3),'b',T,PYR_upt2(:,4),'m',T,PYR_upt2(:,5),'c','LineWidth',2,'MarkerSize',mk_size)
 box off 
text(7,10,'[PYR]=5mM','FontSize',text_size);
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
set(gca,'XTick',(0:2:10),'YTick',(0:30:150));
% hold on
legend('[MAL]=0 mM','[MAL]=0.5 mM','[MAL]=1 mM','[MAL]=5 mM','[MAL]=10 mM','Location','NorthWest')
legend boxoff
ylabel('Mass (nmoles/mg)')
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',text_size)
title('PYR Uptake')
subplot(1,3,2)
plot(Fx1,Fy1,'ro',Fx2,Fy2,'g*',Fx3,Fy3,'b>',Fx4,Fy4,'m+',Fx5,Fy5,'c^',T,MAL_upt2(:,1),'r',T,MAL_upt2(:,2),'g',T,MAL_upt2(:,3),'b',T,MAL_upt2(:,4),'m',T,MAL_upt2(:,5),'c','LineWidth',2,'MarkerSize',mk_size)
title('MAL Uptake')
box off 
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
set(gca,'XTick',(0:2:10),'YTick',(0:30:150));
ylim([-2 150])
xlim([0 11])
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',text_size)

subplot(1,3,3)
plot(Dx1,Dy1,'ro',Dx2,Dy2,'g*',Dx3,Dy3,'b>',Dx4,Dy4,'m+',Dx5,Dy5,'c^',T,CIT_rel2(:,1),'r',T,CIT_rel2(:,2),'g',T,CIT_rel2(:,3),'b',T,CIT_rel2(:,4),'m',T,CIT_rel2(:,5),'c','LineWidth',2,'MarkerSize',mk_size)
box off  
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
set(gca,'XTick',(0:2:10),'YTick',(0:30:150));
hold on
xlim([0 11])
title('CIT Production')
%ylabel('CIT production (nmoles/mg)')
%legend('MAL=0.5 mM','MAL=1 mM','MAL=5 mM','MAL=10 mM')
legend boxoff
xlabel('Time (min)','FontName','Times New Roman', 'FontSize',text_size)
%--------------------------
