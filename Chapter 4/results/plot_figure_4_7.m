% Plot figure 4.7, the figure adjustement depends on the screen
% Run save_sim_M3.m, save_sim_M4.m, save_waic_M3.m and save_waic_M4.m before 
% the variable yp inside simM4.mat may prevent from loading the data
% this variable is not needed for the plot, you can load the other variable
% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig
clear all
close all
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');

S3                      = load('simM3_2.mat');
S4                      = load('simM4_2.mat'); 
W3                      = load('waicM3_2.mat');
W4                      = load('waicM4_2.mat');
X3                      = S3.X;
Y3                      = S3.Y;
X4                      = S4.X;
Y4                      = S4.Y;
ypt                     = S4.ypt;
LL3p                    = W3.LL3p;
LL4p                    = W4.LL4p;

cutp                    = 5001:6441;
Tp                      = length(cutp);
startDate               = datenum('13-05-2014 17:30:00', 'dd-mm-yy HH:MM:SS');
endDate                 = datenum('23-05-2014 17:30:00', 'dd-mm-yy HH:MM:SS'); 
xData                   = linspace(startDate,endDate,Tp-1);

lw                      = 2;
plotheight              = 15;
plotwidth               = 15;
fs                      = 14; % fontSize
fw                      = 'demi'; % fontWeight
fn                      = 'Times New Roman'; % fontName
dcolors                 = linspecer(9,'qualitative');

subplotsx               = 2;
subplotsy               = 1;
leftedge                = 0.7;
rightedge               = 0.2;
topedge                 = 0.3;
bottomedge              = 1;
spacex                  = 0.8;
spacey                  = 0.9;
sub_pos                 = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.65 0.5]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
%----------------------------------------------------------------------
axes('position',sub_pos{1,1});
h4                      = fill(X4,Y4,dcolors(6,:),'edgecolor','none');
set(h4,'facealpha',1)
hold on
h3                      = fill(X3,Y3,dcolors(2,:),'edgecolor','none');
set(h3,'facealpha',1)
hold on
plot(xData,ypt(1:Tp-1),'k','linewidth',1.5)
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
set(gcf,'color','w');
ylabel('°C','fontsize',fs,'FontWeight',fw,'FontName',fn)
% text(1.02,0.02,'May','FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
%----------------------------------------------------------------------
axes('position',sub_pos{2,1});
histogram(LL3p,'FaceColor',dcolors(2,:),'EdgeColor',dcolors(2,:),'EdgeAlpha',0,'Normalization','pdf');
hold on
histogram(LL4p,'FaceColor',dcolors(6,:),'EdgeColor',dcolors(4,:),'EdgeAlpha',0,'Normalization','pdf');
hold on
[f3,xi3]  = ksdensity(LL3p);
plot(xi3,f3,'Color',dcolors(2,:),'Linewidth',lw)
hold on
[f4,xi4]  = ksdensity(LL4p);
plot(xi4,f4,'Color',dcolors(6,:),'Linewidth',lw)
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
