
clear all

clc
format long

%----
text_size=18;
mk_size=12;
line_width=1.5;
%----------------------
%%  Parameter Setup
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
%Test:include R123 as state variables
iR123e=41; iR123m=42; iReBe=43; iRmBm=44;
F_con   =  0.096484;    % kJ mol^{-1} mV^{-1}                          % Faraday 's constant [coulomb/mole]
 %Tem=310.15; %K      37 oC
  Tem=303.15; %K      30 oC
 % Tem=298.15; %K      25 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
 Ve  =   1.9/0.4*1e-3; %L    
Vmito=0.001e-3;% total volume  L (for 1mg mitochondria) 
Vm  = 1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
ROTi=1;
AAi=1;
closed_system=1;% closed system
% % Temperature correction
%activation_energy
Q10=2.5*ones(25,1);
Tem_standard=303.15;  %K  30
Q10_factor=Q10.^((Tem-Tem_standard)/10);

Para=ones(27,1);
 Para(1:25)=1*Q10_factor.*Para(1:25);
 
 substrates=6;  %choose substrate 
  %%  Define t_step and t_final
t_step      =   0.01;   %min
%% Run Simulation
 dTEM=[(303.15-7):2:(303.15+7)];

 O_s3=zeros(5,1);
 O_s2=zeros(5,1);
for i=1:1:length(dTEM)
 ADP_add=30e-6;

 Tem=dTEM(i);

Para=ones(27,1);
 Q10_factor=Q10.^((Tem-Tem_standard)/10);
  Para(1:25)=1*Q10_factor.*Para(1:25);
 dPara=Para;
% Tem=310.15;
time1=2; %Length of State 2  (time before adding ADP)
time2=14;
time3=14;

IC=Set_Initial_Concentrations(substrates,bufferpH);
IC(iPie)=4e-3;
options = odeset('RelTol',1e-10, 'AbsTol',1e-10, 'NormControl','on', ...
          'MaxStep',t_step/5, 'InitialStep',t_step/10, 'MaxOrder',5, ...
          'BDF','on','NonNegative',[1:42]);

[T1,C1] = ode15s(@odeq,[0:t_step:time1],IC,options,substrates,dPara,2);
IC2=C1(end,:); %Intial concentration of State III
IC2(iADPe)=IC2(iADPe)+ ADP_add; %add ADP, Unit(Molar) 

[T2,C2]= ode15s(@odeq,[0:t_step:time2],IC2,options,substrates,dPara,2);


T=[T1; T2(2:end)+time1;];

C=[C1; C2(2:end,:);];

index_O2end1=find(C2(:,iADPe)<0.05*ADP_add,1,'first') %State 3 ends
% 1e9*(Ve+Vm+Vi)*(C1(5,40)-C1(10,40))/(T1(10)-T1(5))%nmols

Tfluxes=zeros(length(C(:,1)),12)';
Rfluxes=zeros(length(C(:,1)),15)';
for istep=1:1:(length(C(1:end,1)))
    RTfluxes(:,istep)=fluxes(C(istep,:),Para);
end
Rfluxes=RTfluxes(1:15,:);
Tfluxes=RTfluxes(16:27,:);
CIV(:,i)=1e9*Rfluxes(14,:);
 O_s2(i)=1e9*(Ve+Vm+Vi)*(C1(time1/t_step/2,iO2)-C1(time1/t_step,iO2))/(T1(time1/t_step)-T1(time1/t_step/2));%nmols
 O_s3(i)=1e9*(Ve+Vm+Vi)*(C2(1,iO2)-C2(index_O2end1,iO2))/(T2(index_O2end1)-T2(1));
Oxy(i,:)=C(:,iO2);
Dfi(i,:)=C(:,idPsi);
NADH(i,:)=C(:,iNADHm);
clear C
end
% Tfluxes=zeros(length(C(:,1)),11)';
% Rfluxes=zeros(length(C(:,1)),14)';
% for istep=1:1:length(C(:,1))
% [Tfluxes(:,istep),Rfluxes(:,istep)]=calculate_flux(C(istep,:),Para);
% end

%%----------------------


set(figure(1),'Units','inches','Position',[0.2 0.1 5 4]) 
set(figure(2),'Units','inches','Position',[5.2 0.1 5 4]) 
set(figure(3),'Units','inches','Position',[0.2 4.1 5 4]) 
set(figure(4),'Units','inches','Position',[5.2 4.1 5 4]) 
set(figure(5),'Units','inches','Position',[10.2 4.1 5 4]) 
set(figure(6),'Units','inches','Position',[10.2 4.1 5 4]) 
%------------------------------------

%------------------------------------------
figure(1)
plot(T,Dfi(1,:),'r',T,Dfi(2,:),'g',T,Dfi(3,:),'b',T,Dfi(4,:),'m',T,Dfi(5,:),'c','LineWidth',line_width)
xlabel('Time (min)')
ylabel('Membrane Potential (mV)')
title('Membrane Potential')
xlim([0 6])
ylim([130 150])
box off
legend(['25' char(176) 'C'],['28' char(176) 'C'],['31' char(176) 'C'],['34' char(176) 'C'],['37' char(176) 'C'])
legend boxoff
  set(gcf,'color','w')
set(gca,'Fontsize',text_size)
figure(2)
plot(T,1e3*NADH(1,:),'r',T,1e3*NADH(2,:),'g',T,1e3*NADH(3,:),'b',T,1e3*NADH(4,:),'m',T,1e3*NADH(5,:),'c','LineWidth',line_width)
box off
xlabel('Time (min)')
ylabel('Concentration (mM)')
xlim([0 6])
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
legend(['25' char(176) 'C'],['28' char(176) 'C'],['31' char(176) 'C'],['34' char(176) 'C'],['37' char(176) 'C'])
legend boxoff
title('NADH','FontSize',text_size)

figure(3)
plot(T,Oxy(1,:),'r',T,Oxy(2,:),'g',T,Oxy(3,:),'b',T,Oxy(4,:),'m',T,Oxy(5,:),'c','LineWidth',line_width)
xlabel('Time (min)')
ylabel('Concentration (M)')
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
legend(['25' char(176) 'C'],['28' char(176) 'C'],['31' char(176) 'C'],['34' char(176) 'C'],['37' char(176) 'C'])
legend boxoff
title('Oxygen Concentration','FontSize',text_size)
Tem_index=[23 28 30];
State2_PM=[6.35/2 7.1 12.6];
State3_PM=[35.46/2 42.8 35];
S2_SE_PM=[0.08 0 1.7 ]
S3_SE_PM=[2.24 0  2.4]
% State2_SUC=[10.61/2 53/2.1 26.3 17];  % SUC 28 degree Fisher et al.
% State3_SUC=[60.34/2 53 68.9 55];
State2_SUC=[10.61/2 53/2.1 26.3];  % SUC 28 degree Fisher et al.
State3_SUC=[60.34/2 53 68.9];
S2_SE_SUC=[0.9 9/2.1 4 ]
S3_SE_SUC=[3.83 9 4]
figure(4)
if substrates==6
plot(dTEM-273.15,O_s2,'b',dTEM-273.15,O_s3,'r',Tem_index,State2_PM,'b*',Tem_index,State3_PM,'ro','LineWidth',line_width,'MarkerSize',8)
hold on 
h2=errorbar(Tem_index,State2_PM,S2_SE_PM,'d','marker','square','markersize',8,...
                 'markeredgecolor','b',...
                'color','b','linewidth',1.5)
            hold on
 h2=errorbar(Tem_index,State3_PM,S3_SE_PM,'d','marker','square','markersize',8,...
                 'markeredgecolor','r',...
                'color','r','linewidth',1.5)
            hold off
                 title('OCR (PYR+MAL)','FontSize',text_size)
elseif substrates==3
 plot(dTEM-273.15,O_s2,'b',dTEM-273.15,O_s3,'r',Tem_index,State2_SUC,'b*',Tem_index,State3_SUC,'ro','LineWidth',line_width,'MarkerSize',8)
hold on
 h2=errorbar(Tem_index,State2_SUC,S2_SE_SUC,'d','marker','square','markersize',8,...
                 'markeredgecolor','b',...
                'color','b','linewidth',1.5)
       
hold on
 h2=errorbar(Tem_index,State3_SUC,S3_SE_SUC,'d','marker','square','markersize',8,...
                 'markeredgecolor','r',...
                'color','r','linewidth',1.5)
            hold off
                 title('OCR (SUC)','FontSize',text_size)
end
xlim([23 37])
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
legend('State 2','State 3')
legend boxoff
xlabel(['Temperature ' char(176) 'C'],'FontName','Times New Roman', 'FontSize',text_size)
ylabel('Respiration (nmol/min/mg)')
% title('OCR','FontSize',text_size)
figure(5)
plot(dTEM-273.15,O_s3./O_s2,'b','LineWidth',line_width)
xlabel(['Temperature ' char(176) 'C'],'FontName','Times New Roman', 'FontSize',text_size)
title('RCI','FontSize',text_size)
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
%-------------------------------------

figure(6)
plot(T,CIV(:,1),'r',T,CIV(:,2),'g',T,CIV(:,3),'b',T,CIV(:,4),'c',T,CIV(:,5),'m','LineWidth',line_width)
set(gcf,'color','w')
set(gca,'Fontsize',text_size)
title('OCR','FontSize',text_size)
legend(['25' char(176) 'C'],['28' char(176) 'C'],['31' char(176) 'C'],['34' char(176) 'C'],['37' char(176) 'C'])
legend boxoff

xlabel('Time (min)','FontName','Times New Roman', 'FontSize',text_size)
ylabel('OCR (nmol/min/mg)')


% plot(dTEM-273.15,State2_PM,'o',Tem,State3_PM,'*')
% ylim([0 55])
% title('PYR+MAL')
%  figure(3)
% plot(Tem,State2_SUC,'o',Tem,State3_SUC,'*')
% ylim([0 70])
% title('SUC')
% figure(1)