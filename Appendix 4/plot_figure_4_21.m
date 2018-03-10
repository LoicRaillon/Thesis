% Generate Figure 4.21, pvalue of t-test
% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig

clear all
close all
clc
addpath(userpath, 'utilities');

plotheight      = 15;
plotwidth       = 15;
subplotsx       = 2;
subplotsy       = 1;
leftedge        = 0.9;
rightedge       = 0.5;
topedge         = 0.6;
bottomedge      = 1.4;
spacex          = 1;
spacey          = 1;
sub_pos         = subplot_pos(plotwidth,plotheight,leftedge,rightedge,...
                    bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
lw              = 1.5;
plotheight      = 15;
plotwidth       = 15;
fs              = 14;                                                       % fontSize
fw              = 'demi';                                                   % fontWeight
fn              = 'Times New Roman';                                        % fontName
figure('visible','on','units','normalized', ...
    'outerposition',[0.5 0.04 0.5 0.45]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');

%% -------------------------------------------------------------------------
% Left Plot 
% -------------------------------------------------------------------------
N               = 1e6;                                                      % Number of sample for tpdf plot
dof             = 2111-19;                                                  % degree of freedom: Nsamples - Nparameters
alpha           = 0.05;                                                     % significance level 
ra              = tinv(1-alpha/2,dof);                                      % left critical region t-student pdf
la              = tinv(alpha/2,dof);                                        % right critical region t-student pdf 
mle             = 2.37e-2;                                                  % maximum likelihood estimate
se              = 5.45e-3;                                                  % standard deviation of maximum likelihood estimate
z               = mle/se;                                                   % t-value (4.54)
pval1           = 2*(1-tcdf(abs(z),dof));                                   % pvalue (4.56)

x               = sort([linspace(-5,5,N) ra la z]);                         % x grid
y0              = tpdf(x,dof);                                              % t-student pdf
idx             = find(x >= la & x <= ra);                                  % find critical region (black)
yc              = ones(1,length(y0));                                       % build mask
yc(idx)         = 0;                                                        % ...
idx2            = find(x >= -z & x <= z);                                   % find pvalue region (red)
yz              = ones(1,length(y0));                                       % build mask
yz(idx2)        = 0;                                                        % ...

axes('position',sub_pos{1,1});
AR              = area(x,yc.*y0); hold on;                                  % fill critical region with black
AR.EdgeColor    = 'k';
AR.EdgeAlpha    = 1;
AR.FaceColor    = 'k';
AR.FaceAlpha    = 1;
AR2             = area(x,yz.*y0); hold on;                                  % fill pvalue region with red
AR2.EdgeColor   = 'k';
AR2.EdgeAlpha   = 1;
AR2.FaceColor   = 'r';
AR2.FaceAlpha   = 1;
plot(x,y0,'k','Linewidth',1); hold on;
xx              = get(gca,'Xtick');
plot([xx(1) -z],[3e-3 3e-3],'r','Linewidth',2); hold on;
plot([z xx(end)],[3e-3 3e-3],'r','Linewidth',2); hold on;
axis tight
set(gca,'Xtick',[-4 -2 0 2 4]);
set(gca,'Ytick',[0 0.1 0.2 0.3 0.4]);
set(gcf,'color','w');
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fs)
set(gca,'FontName',fn)
set(gca,'Linewidth',lw)

%% -------------------------------------------------------------------------
% Right Plot, similar to above except that the std is doubled 
% -------------------------------------------------------------------------
z               = mle/(2*se); 
pval2           = 2*(1-tcdf(abs(z),dof));

x               = sort([linspace(-4,4,N) ra la z]);
y0              = tpdf(x,dof);
idx             = find(x >= la & x <= ra); 
yc              = ones(1,length(y0));
yc(idx)         = 0; 
idx2            = find(x >= -z & x <= z); 
yz              = ones(1,length(y0));
yz(idx2)        = 0; 

axes('position',sub_pos{2,1});
AR              = area(x,yc.*y0); hold on;
AR.EdgeColor    = 'k';
AR.EdgeAlpha    = 1;
AR.FaceColor    = 'k';
AR.FaceAlpha    = 1;
AR2             = area(x,yz.*y0); hold on;
AR2.EdgeColor   = 'k';
AR2.EdgeAlpha   = 1;
AR2.FaceColor   = 'r';
AR2.FaceAlpha   = 1;
plot(x,y0,'k','Linewidth',1); hold on;
axis tight
set(gca,'Ytick',[0 0.1 0.2 0.3 0.4]);
set(gcf,'color','w');
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fs)
set(gca,'FontName',fn)
set(gca,'Linewidth',lw)
