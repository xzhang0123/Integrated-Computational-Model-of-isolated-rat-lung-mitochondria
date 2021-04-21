%plot 3D surface figure with saved data
close all
clear all
clc
text_size=8;
text_size2=14;
title_size=15;
line_width=2;
load data0418.mat
%-------------PYR-----------
set(figure(1),'Units','inches','Position',[2 1 5*3 4*3]) 
set(figure(2),'Units','inches','Position',[1 1 1.1*5*4 2.2*3.5]) 
set(figure(3),'Units','inches','Position',[1 1 1.1*5*4 2.2*3.5*1]) 
set(figure(4),'Units','inches','Position',[1 1 1.1*5*4 2.2*3.5*1]) 
figure(1)
subplot(3,3,1)
surf(1e3*MALstep,1e3*PYRstep,PYRm,'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10]) 
zlim([0 11])
text(0.5,9,6,'PYR','Fontsize',text_size2)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Concentration (mM)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
%title('PYRm and MALm','Fontsize',title_size)
%-------------------------------
subplot(3,3,2)
% surf(1e3*MALstep,1e3*PYRstep,ACOAm./(ACOAm+COAm));
surf(1e3*MALstep,1e3*PYRstep,ACOAm,'EdgeColor','none','FaceAlpha',.8);colormap jet
hold on
surf(1e3*MALstep,1e3*PYRstep,COAm,'EdgeColor','none','FaceAlpha',.8);colormap jet;
hold off
xlim([0 10])
ylim([0 10])
zlim([0 1.5])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
text(0,9,0.3,'COA','Fontsize',text_size2)
text(8,2,0.3,'ACOA','Fontsize',text_size2)

 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Concentration (mM)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
subplot(3,3,3)
surf(1e3*MALstep,1e3*PYRstep,CITm,'EdgeColor','none','FaceAlpha',.8);colormap jet;
text(0,9,8,'CIT','Fontsize',text_size2)
xlim([0 10])
ylim([0 10])
zlim([0 15])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Concentration (mM)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
%-----------------------------------------------------

subplot(3,3,4)
surf(1e3*MALstep,1e3*PYRstep,SUCm,'EdgeColor','none','FaceAlpha',.8)
hold on
surf(1e3*MALstep,1e3*PYRstep,FUMm,'EdgeColor','none','FaceAlpha',.8)
hold on
surf(1e3*MALstep,1e3*PYRstep,MALm,'EdgeColor','none','FaceAlpha',.8);colormap jet;
hold off
text(0.3,10,5,'FUM','Fontsize',text_size2)
text(0.3,10,1,'SUC','Fontsize',text_size2)
text(0.3,10,10,'MAL','Fontsize',text_size2)
xlim([0 10])
ylim([0 10])
zlim([0 15])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Concentration (mM)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
%-----------------------------------------------------
subplot(3,3,5)
surf(1e3*MALstep,1e3*PYRstep,NADHm./(NADHm+NADm),'EdgeColor','none','FaceAlpha',.8);colormap jet;
hold on
surf(1e3*MALstep,1e3*PYRstep,FADH2m./(FADH2m+FADm),'EdgeColor','none','FaceAlpha',.8);colormap jet;
hold off
xlim([0 10])
ylim([0 10])
text(0.5,9,1,'NADH','Fontsize',text_size2)
text(0.5,9,0.3,'FADH_2','Fontsize',text_size2)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Ratio','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
subplot(3,3,6)
surf(1e3*MALstep,1e3*PYRstep,UQH2m./(UQm+UQH2m),'EdgeColor','none','FaceAlpha',.8);colormap jet;
hold on
xlim([0 10])
ylim([0 10])
text(0.5,10,8e-3,'UQH_2','Fontsize',text_size2)
%surf(1e3*MALstep,1e3*PYRstep,UQm,'FaceColor','r','EdgeColor','none','FaceAlpha',.8)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Ratio','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
%-----------------------------------------------------
subplot(3,3,7)
surf(1e3*MALstep,1e3*PYRstep,CytCredm./(CytCredm+CytCoxim),'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
%hold on
%surf(1e3*MALstep,1e3*PYRstep,CytCoxim,'FaceColor','r','EdgeColor','none','FaceAlpha',.8)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
text(0.5,10,0.22,'CytC(ratio)','Fontsize',text_size2)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Ratio','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
%-----------------------------------------------------
subplot(3,3,8)
surf(1e3*MALstep,1e3*PYRstep,ATPm./(ATPm+ADPm),'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
text(0.5,9,0.09,'ATP','Fontsize',text_size2)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Ratio','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
%-----------------------------------------------------
subplot(3,3,9)
surf(1e3*MALstep,1e3*PYRstep,dPsi,'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
zlim([145 150])
text(0.5,10,145,'d\Psi','Fontsize',text_size2)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','left','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Membrane Potential (mV)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
% %% 
% %----------------------------------------------------------
%  text_size=9;
figure(2)
subplot(2,4,1)
surf(1e3*MALstep,1e3*PYRstep,PYRH,'EdgeColor','none','FaceAlpha',.8);colormap jet;
text(0.5,9,18,'PYRH','Fontsize',text_size2,'FontWeight','bold')
xlim([0 10])
ylim([0 10])
zlim([2 15])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-------------------------------
subplot(2,4,2)
surf(1e3*MALstep,1e3*PYRstep,TCC,'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
zlim([2 15])
text(0.5,9,18,'TCC','Fontsize',text_size2,'FontWeight','bold')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%--------------------------------
subplot(2,4,3)
surf(1e3*MALstep,1e3*PYRstep,DCCMAL,'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
zlim([-2 3])
text(0.5,9,18/13*3,'DCC(MAL)','Fontsize',text_size2,'FontWeight','bold')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
subplot(2,4,4)
surf(1e3*MALstep(5:end),1e3*PYRstep,DCCSUC(:,5:end),'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
zlim([-1 0.5])
text(0.5,9,18/13*0.5,'DCC(SUC)','Fontsize',text_size2,'FontWeight','bold')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');



figure(3)
subplot(2,4,1)
surf(1e3*MALstep,1e3*PYRstep,PYRDH,'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
zlim([0 15])
text(0.5,9,17,'PDH','Fontsize',text_size2,'FontWeight','bold')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
subplot(2,4,2)
surf(1e3*MALstep,1e3*PYRstep,MALDH,'EdgeColor','none','FaceAlpha',.8);colormap jet;
text(0.5,9,17,'MDH','Fontsize',text_size2,'FontWeight','bold')
xlim([0 10])
ylim([0 10])
zlim([0 15])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');


%-----------------------------------------------------
subplot(2,4,3)
surf(1e3*MALstep,1e3*PYRstep,CITD,'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
zlim([0 15])
text(0.5,9,17,'CITDH','Fontsize',text_size2,'FontWeight','bold')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
subplot(2,4,4)
surf(1e3*MALstep,1e3*PYRstep,AKGDH,'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
zlim([0 15])
text(0.5,9,17,'AKGDH','Fontsize',text_size2,'FontWeight','bold')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------

figure(4)


subplot(2,4,1)
surf(1e3*MALstep,1e3*PYRstep,CI,'EdgeColor','none','FaceAlpha',.8);colormap jet;
xlim([0 10])
ylim([0 10])
zlim([15 30])
% title('CI')
text(0.5,9,22,'CI','Fontsize',text_size2)
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
sp=5;
subplot(2,4,2)
surf(1e3*MALstep(sp:end),1e3*PYRstep,CIII(:,sp:end),'EdgeColor','none','FaceAlpha',.8);colormap jet;
text(0.5,9,23,'CIII','Fontsize',text_size2)
xlim([0 10])
ylim([0 10])
zlim([28 30])
% title('CIII')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
subplot(2,4,3)
surf(1e3*MALstep(sp:end),1e3*PYRstep,CIV(:,sp:end),'EdgeColor','none','FaceAlpha',.8);colormap jet;
text(0.5,9,23,'CIV','Fontsize',text_size2)
xlim([0 10])
ylim([0 10])
zlim([28 30])
% title('CIV')
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
subplot(2,4,4)
surf(1e3*MALstep(sp:end),1e3*PYRstep,leak(:,sp:end),'EdgeColor','none','FaceAlpha',.8);colormap jet;
text(0.5,9,225,'Leak','Fontsize',text_size2)
% title('Leak')
xlim([0 10])
ylim([0 10])
zlim([280 300])
set(gcf,'color','w')
set(gca,'Fontsize',text_size,'LineWidth',line_width)
 xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
 ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

%-----------------------------------------------------
% 
% figure(5)
% 
% subplot(1,4,1)
% 
% 
% surf(1e3*MALstep,1e3*PYRstep,CII,'EdgeColor','none','FaceAlpha',.8);colormap jet;
% % zlim([19 22])
% title('CII')
% text(0.5,9,9,'CI','Fontsize',text_size2)
% set(gcf,'color','w')
% set(gca,'Fontsize',text_size,'LineWidth',line_width)
%  xlabel('[PYR] (mM)','Rotation',15,'HorizontalAlignment','center','VerticalAlignment','middle'); 
%  ylabel('[MAL] (mM)','Rotation',330,'HorizontalAlignment','center','VerticalAlignment','middle'); 
% zlabel('Flux (nmol/min/mg)','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
% 
% %-----------------------------------------------------