function dxdt = odeq(t,IC,substrates,Para1,state)

%   The purpose of this function is to establish the rate equations for
%   chemical reactions of molecular compounds and transport of these
%   compounds between compartments in the lung cell metabolism model.
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
 % Tem=310.15; %K  

buffer_m= 2.303*10^(-7.6)/43; 

%buffer capacity =0.0043 M in matrix
%----------------------------------------------
rtfluxes=fluxes(IC,Para1);
J_PDH=rtfluxes(1);
J_CITS=rtfluxes(2);
J_CITD=rtfluxes(3);
J_AKGD=rtfluxes(4);
J_SCAS=rtfluxes(5);
J_NDK=rtfluxes(6);
J_SUCD=rtfluxes(7);
J_FUM=rtfluxes(8);
J_MALD=rtfluxes(9);
J_GOT=rtfluxes(10);
J_CI=rtfluxes(11);
J_CII=rtfluxes(12);
J_CIII=rtfluxes(13);
J_CIV=rtfluxes(14);
J_CV=rtfluxes(15);
T_SUC_Pi=rtfluxes(16);
T_MAL_Pi=rtfluxes(17);
T_MAL_aKG=rtfluxes(18);
T_MAL_HCIT=rtfluxes(19);
T_PYRH=rtfluxes(20);
T_PIC=rtfluxes(21);
T_ANT=rtfluxes(22);
T_GLUH=rtfluxes(23);
T_ASP_GLU=rtfluxes(24);
T_HLEAK=rtfluxes(25);
T_FUM_Pi=rtfluxes(26);
T_R123=rtfluxes(27);
%--------------------
if state==3&&t>4
    T_HLEAK=10*T_HLEAK;
end
   
% J_CIII=0.5*J_CIII;
% J_CIV=0.5*J_CIV;
R123e=IC(41);
R123m=IC(42);
ReBe=IC(43);
RmBm=IC(44);

%% Net Equations

dxdt(iPie)= ( +T_SUC_Pi-T_PIC+T_MAL_Pi)/Ve;
dxdt(iADPe) = (-T_ANT)/Ve;
dxdt(iATPe) = (+T_ANT)/Ve;
dxdt(iPYRe) = (-T_PYRH)/Ve;
dxdt(iMALe) = (-T_MAL_Pi-T_MAL_aKG-T_MAL_HCIT)/Ve;
dxdt(iCITe) = (+T_MAL_HCIT)/Ve;
dxdt(iaKGe) = (+T_MAL_aKG)/Ve;
dxdt(iSUCe) = (-T_SUC_Pi)/Ve;
dxdt(iFUMe) = (T_FUM_Pi)/Ve;
dxdt(iGLUe) = (-T_GLUH+T_ASP_GLU)/Ve;
dxdt(iASPe) = (-T_ASP_GLU)/Ve;
%Mitochondria
dxdt(iGLUm) = (J_GOT - T_ASP_GLU + T_GLUH )/Vm;
dxdt(iASPm) = (-J_GOT+T_ASP_GLU)/Vm;
dxdt(iPYRm) = (T_PYRH - J_PDH)/Vm;
dxdt(iOXAm)= (-J_CITS + J_MALD+J_GOT)/Vm;
dxdt(iCITm) = (-T_MAL_HCIT + J_CITS - J_CITD)/Vm;
 %dxdt(iCITm)=0;
dxdt(iaKGm)  = (J_CITD - J_GOT - J_AKGD - T_MAL_aKG)/Vm;
dxdt(iSCAm)  = (J_AKGD - J_SCAS)/Vm;
dxdt(iSUCm)  = (+T_SUC_Pi + J_SCAS - J_SUCD)/Vm;
dxdt(iFUMm)  =(-T_FUM_Pi+J_SUCD-J_FUM)/Vm;
dxdt(iMALm)  = (J_FUM - J_MALD  +T_MAL_Pi + T_MAL_aKG + T_MAL_HCIT)/Vm;
dxdt(iNADm) = (-J_PDH - J_CITD - J_AKGD - J_MALD + J_CI)/Vm;
dxdt(iNADHm) = -dxdt(iNADm);
dxdt(iUQm)=(-J_CI-J_CII+J_CIII)/Vm;
dxdt(iUQH2m)=-dxdt(iUQm);
dxdt(iCytCoxi)=(-2*J_CIII+2*J_CIV)/Vi;
dxdt(iCytCred)=-dxdt(iCytCoxi);
dxdt(iADPm) = (T_ANT - J_NDK - J_CV )/Vm;
dxdt(iATPm) = (J_CV + J_NDK - T_ANT )/Vm;
dxdt(iGDPm) = (-J_SCAS+J_NDK)/Vm;
dxdt(iGTPm) = -dxdt(iGDPm);
dxdt(iCOAm) = (- J_PDH+J_CITS + J_SCAS - J_AKGD )/Vm;
 dxdt(iACOAm)= (J_PDH-J_CITS)/Vm;

dxdt(iPim)  = (-T_SUC_Pi - T_MAL_Pi - J_SCAS - J_CV +T_PIC  )/Vm;
dxdt(iFADm)  = (J_CII-J_SUCD)/Vm;
dxdt(iFADH2m)= -dxdt(iFADm) ;
%beta=beta'/([H+]*2.303)
%buffering capacity in cytosol for protons=6.65 mmol/pH
%buffering capacity in mitochondria       =25 mmol/pH
dxdt(iHm)   = (buffer_m*( -T_MAL_HCIT+T_PYRH+ T_PIC+T_GLUH-T_ASP_GLU+T_HLEAK- J_PDH+2*J_CITS +(2-2)* J_CITD+(1-2)*J_AKGD+J_SCAS+ J_MALD - (4+1)*J_CI - 2*J_CIII-4*J_CIV+(3-1)*J_CV))/Vm;
% dxdt(iHe)   = buffer_c*(4*J_CI+4*J_CIII+2*J_CIV-3*J_CV-T_HLEAK+ T_ASP_GLU+T_MAL_HCIT-T_PYRH-T_PIC-T_GLUH+T_ASP_GLU)/Ve;
dxdt(iHe)   =0;
%6.75e-3 mmol/mV is the capacitance of the IMM, or 6.75e3 umol/V
Cimm=1*6.75*1e-6*(Vm+Vi);       
%6.75*1e-6mol/(L mito*mV). Mito volume is 1.3e-6 L in Jason's paper
dxdt(idPsi)=1/Cimm*(4*J_CI+2*J_CIII+4*J_CIV-3*J_CV-T_HLEAK-T_ANT-T_R123-T_ASP_GLU);
dxdt(iO2)=closed_system*0.5*(-J_CIV)/(Vm+Ve+Vi);
%--------------------
%Test:include R123 as state variables
k1p=0*4.3e-3;
k1n=0*8.9e-2;
k2p=0*4.3e-2;
k2n=0*4.45e-2;
% dxdt(iR123e)=0;
% dxdt(iR123m)=0;
dxdt(iRmBm)=k2p*R123m-k2n*RmBm;
dxdt(iReBe)=k1p*R123e-k1n*ReBe;
% %----------------------
 dxdt(iR123m)=T_R123/(Vm*2.255)-k2p*R123m+k2n*RmBm;
 dxdt(iR123e)=-T_R123/Ve-k1p*R123e+k1n*ReBe;
% isolated  ETC
% dxdt(iPYRe)=0;
% dxdt(iMALe)=0;
% dxdt(iCITe)=0;

%--------------
dxdt=dxdt';
%% Output Matrix

end

