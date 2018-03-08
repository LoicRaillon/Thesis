% Plot figure 4.7, the figure adjustement depends on the screen
% Run Profile_ML_MAP_M3 before
% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig

clear all
close all
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');

lw              = 2;
plotheight      = 15;
plotwidth       = 15;
fs              = 14; % fontSize
fw              = 'demi'; % fontWeight
fn              = 'Times New Roman'; % fontName

subplotsx       = 2;
subplotsy       = 1;
leftedge        = 0.7;
rightedge       = 0.8;
topedge         = 0.5;
bottomedge      = 1;
spacex          = 1.4;
spacey          = 0.9;
sub_pos         = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

figure('visible','on','units','normalized','outerposition',[0 0.0370 0.5000 0.9630/2]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');

imp             = load('PL_ML_MAP_sw33_bis.mat');
g               = imp.g;
pML             = imp.pML;
pMAP            = imp.pMAP;

grid            = linspace(min(g),max(g),1e2)';
pML2            = pchip(g,pML,grid);
pMAP2           = pchip(g,pMAP,grid);
mpML2           = max(pML2);
mpMAP2          = max(pMAP2);
idxML           = find(pML2 == mpML2);
idxMAP          = find(pMAP2 == mpMAP2);

axes('position',sub_pos{1,1});
plot(grid,pML2 - mpML2,'b',grid,pMAP2 - mpMAP2,'r','Linewidth',lw); hold on;
plot(grid(idxML),0,'ko','Linewidth',lw,'MarkerFaceColor','g'); hold on;
plot(grid(idxMAP),0,'ko','Linewidth',lw,'MarkerFaceColor','g'); hold on;
plot([grid(1) grid(end)],[-1.92 -1.92],'--k','Linewidth',lw)
ylim([-3 0])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xx              = get(gca,'Xtick');
nx              = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.015,0.03,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
text(0.8,0.95,'\textbf{$\sigma_{w_{33}}$}','FontSize',fs+4,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')

imp             = load('PL_ML_MAP_Ci_bis.mat');
g               = imp.g;
pML             = imp.pML;
pMAP            = imp.pMAP;

grid            = linspace(min(g),max(g),1e2)';
pML2            = pchip(g,pML,grid);
pMAP2           = pchip(g,pMAP,grid);
mpML2           = max(pML2);
mpMAP2          = max(pMAP2);
idxML           = find(pML2 == mpML2);
idxMAP          = find(pMAP2 == mpMAP2);

axes('position',sub_pos{2,1});
plot(grid,pML2 - mpML2,'b',grid,pMAP2 - mpMAP2,'r','Linewidth',lw); hold on;
plot(grid(idxML),0,'ko','Linewidth',lw,'MarkerFaceColor','g'); hold on;
plot(grid(idxMAP),0,'ko','Linewidth',lw,'MarkerFaceColor','g'); hold on;
plot([grid(1) grid(end)],[-1.92 -1.92],'--k','Linewidth',lw)
ylim([-3 0])
xlim([6.57e-3 6.8e-3])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xx              = get(gca,'Xtick');
nx              = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.015,0.03,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
text(0.8,0.95,'\textbf{$C_i$}','FontSize',fs+4,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
