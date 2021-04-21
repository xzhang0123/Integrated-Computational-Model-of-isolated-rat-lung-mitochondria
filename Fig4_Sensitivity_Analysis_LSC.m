clear all
clc
 
global k0 Tem C F_con R_con Ve Vm Vi ROTi AAi closed_system...
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
 Tem=303.15; %K      30 oC
R_con  = 8.314e-3;   %gas constant [kJ/K/mol]
ROTi=1;
AAi=1;
%volumes
Ve  =   1.9/0.4*1e-3; %L  
Vmito=0.65*3.67e-6;% total volume  L (for 1mg mitochondria) 
Vmito=1e-6;% total volume  L (for 1mg mitochondria) 

Vm  = 1e-6;
Vi  =   0.1*Vm;  %L  (10% of mito volume)
%---------------------
closed_system=1;% 
oxi_N8=xlsread('O2_results_n8.xlsx',1);
fs=4;
PM_mean=oxi_N8(44:fs:end,1)/1e6;
SUC_mean=oxi_N8(44:fs:end,3)/1e6;
PM_mean=PM_mean(1:33);
SUC_mean=SUC_mean(1:33);
time_mean=1:3*fs:(3*fs*length(PM_mean));
time_mean=time_mean/60;
%-----------------------------------
oxy_PM=PM_mean;
oxy_SUC=SUC_mean;
size_index=dlmread('text_size.txt');

tic

Ay1=xlsread('TCA1975.xlsx','B3:B7');
Ay2=xlsread('TCA1975.xlsx','B9:B13');
Ay3=xlsread('TCA1975.xlsx','B15:B19');
Ay4=xlsread('TCA1975.xlsx','B21:B25');
Ay5=xlsread('TCA1975.xlsx','B27:B31');

By1=xlsread('TCA1975.xlsx','E3:E7');
By2=xlsread('TCA1975.xlsx','E9:E13');
By3=xlsread('TCA1975.xlsx','E15:E19');
By4=xlsread('TCA1975.xlsx','E21:E25');
By5=xlsread('TCA1975.xlsx','E27:E31');

Cy1=xlsread('TCA1975.xlsx','H3:H7');
Cy2=xlsread('TCA1975.xlsx','H9:H13');
Cy3=xlsread('TCA1975.xlsx','H15:H19');
Cy4=xlsread('TCA1975.xlsx','H21:H25');
Cy5=xlsread('TCA1975.xlsx','H27:H31');

Dy1=xlsread('TCA1975.xlsx','K3:K7');
Dy2=xlsread('TCA1975.xlsx','K9:K13');
Dy3=xlsread('TCA1975.xlsx','K15:K19');
Dy4=xlsread('TCA1975.xlsx','K21:K25');
Dy5=xlsread('TCA1975.xlsx','K27:K31');

Ey1=xlsread('TCA1975.xlsx','M3:M7');
Ey2=xlsread('TCA1975.xlsx','M9:M13');
Ey3=xlsread('TCA1975.xlsx','M15:M19');
Ey4=xlsread('TCA1975.xlsx','M21:M25');
Ey5=xlsread('TCA1975.xlsx','M27:M31');

Fy1=xlsread('TCA1975.xlsx','O3:O7');
Fy2=xlsread('TCA1975.xlsx','O9:O13');
Fy3=xlsread('TCA1975.xlsx','O15:O19');
Fy4=xlsread('TCA1975.xlsx','O21:O25');
Fy5=xlsread('TCA1975.xlsx','O27:O31');
toc
dataTCA=[Ay1;Ay2;Ay3;Ay4;Ay5;By1;By2;By3;By4;By5;Cy1;Cy2;Cy3;Cy4;Cy5;Dy1;Dy2;Dy3;Dy4;Dy5;Ey1;Ey2;Ey3;Ey4;Ey5;Fy1;Fy2;Fy3;Fy4;Fy5] ; %150 by1
%-------------------ydata simulation2-----
t_step=1;
T1=0:t_step:3;
T2=0:t_step:3;
 meanOCR_SUC=[26.3 68.9]; %nmoles/min/mg
 meanOCR_PYRMAL=[12.6 35];
 meanOCR_GM=[22.8 52.15 41.62/2];
 meanOCR_AKGPYR=[10.4 30.1];
 meanOCR_CITPYR=[12.8 42.3];
 initial_o2=296e-6;
    EXP_CITPYR1= initial_o2-meanOCR_CITPYR(1)*1e-9/(Vm+Vi+Ve)*T1;  %Unit is molar
    EXP_CITPYR2=initial_o2-meanOCR_CITPYR(1)*1e-9/(Vm+Vi+Ve)*3-meanOCR_CITPYR(2)*1e-9/(Vm+Vi+Ve)*T2;
EXP_CITPYR=1e6*[EXP_CITPYR1 EXP_CITPYR2(2:end)];
    EXP_SUC1= initial_o2-meanOCR_SUC(1)*1e-9/(Vm+Vi+Ve)*T1;
    EXP_SUC2=initial_o2-meanOCR_SUC(1)*1e-9/(Vm+Vi+Ve)*3-meanOCR_SUC(2)*1e-9/(Vm+Vi+Ve)*T2;
EXP_SUC=1e6*[EXP_SUC1 EXP_SUC2(2:end)];
    EXP_AKGPYR1= initial_o2-meanOCR_AKGPYR(1)*1e-9/(Vm+Vi+Ve)*T1;
    EXP_AKGPYR2=initial_o2-meanOCR_AKGPYR(1)*1e-9/(Vm+Vi+Ve)*3-meanOCR_AKGPYR(2)*1e-9/(Vm+Vi+Ve)*T2;
EXP_AKGPYR=1e6*[EXP_AKGPYR1 EXP_AKGPYR2(2:end)];
    EXP_PYRMAL1= initial_o2-meanOCR_PYRMAL(1)*1e-9/(Vm+Vi+Ve)*T1;
    EXP_PYRMAL2=initial_o2-meanOCR_PYRMAL(1)*1e-9/(Vm+Vi+Ve)*3-meanOCR_PYRMAL(2)*1e-9/(Vm+Vi+Ve)*T2;
EXP_PYRMAL=1e6*[EXP_PYRMAL1 EXP_PYRMAL2(2:end)];
data_OCR=[EXP_CITPYR EXP_SUC EXP_AKGPYR EXP_PYRMAL];
% end of ydata simulation 2 ------------------------------------
Para=ones(27,1);
ydata=[dataTCA./length(dataTCA);1e6*PM_mean(1:33)./33;1e6*SUC_mean(1:33)./33; data_OCR(1:20)'./(length(data_OCR(1:20)))]
for i=1:1:length(Para)

    Para1_loop=Para; Para2_loop=Para;
    Para1_loop(i)=1.01*Para(i);
    Para2_loop(i)=0.99*Para(i);
SSE0(i)=Calculate_SSE(Para,ydata);
SSE1(i)=Calculate_SSE(Para1_loop,ydata);
SSE2(i)=Calculate_SSE(Para2_loop,ydata);
%LSC(i)=(SSE1(i)-SSE2(i))/(0.02*SSE0(i));
end
LSC=abs((SSE1-SSE2)./(0.02*SSE0));
LSC=LSC'