% Generate Figure 4.22, pvalue of Likelihood ratio test
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

%% ------------------------------------------------------------------------
% Left Plot dof = 3
% -------------------------------------------------------------------------
N               = 1e6;                                                      % Number of sample for tpdf plot
dof             = 3;                                                        % Np(M-Msub)
alpha           = 0.05;                                                     % significance level 
ra              = chi2inv(1-alpha,dof);                                     % right critical region t-student pdf
Msub            = 4.2848e+03;                                               % ML large model
M               = 4.2870e+03;                                               % ML sub model
LR              = -2*(Msub - M);                                            % Likelihood ratio (4.59)
pval1           = 1 - pchisq(LR, dof);                                      % pvalue 

x               = sort([linspace(0,15,N) ra LR]);                           % x grid
y0              = chi2pdf(x,dof);                                           % chi2pdf
idx             = find(x <= ra);                                            % find critical region (black)
yc              = ones(1,length(y0));                                       % build mask
yc(idx)         = 0;                                                        % ...
idx2            = find(x <= LR);                                            % find pvalue region (red)
yz              = ones(1,length(y0));                                       % build mask
yz(idx2)        = 0;                                                        % ...

axes('position',sub_pos{1,1});
AR2             = area(x,yz.*y0); hold on;                                  % fill pvalue region with red
AR2.EdgeColor   = 'k';
AR2.EdgeAlpha   = 1;
AR2.FaceColor   = 'r';
AR2.FaceAlpha   = 1;
AR              = area(x,yc.*y0); hold on;                                  % fill critical region with black
AR.EdgeColor    = 'k';
AR.EdgeAlpha    = 1;
AR.FaceColor    = 'k';
AR.FaceAlpha    = 1;
plot(x,y0,'k','Linewidth',1); hold on;
axis tight
set(gcf,'color','w');
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fs)
set(gca,'FontName',fn)
set(gca,'Linewidth',lw)

%% ------------------------------------------------------------------------
% Right Plot dof = 2, similar to above, except the dof
% -------------------------------------------------------------------------
dof             = 2; 
ra              = chi2inv(1-alpha,dof); 
pval2           = 1 - pchisq(LR, dof);

x               = sort([linspace(0,15,N) ra LR]);
y0              = chi2pdf(x,dof);
idx             = find(x <= ra); 
yc              = ones(1,length(y0));
yc(idx)         = 0; 
idx2            = find(x <= LR); 
yz              = ones(1,length(y0));
yz(idx2)        = 0; 

axes('position',sub_pos{2,1});
AR2             = area(x,yz.*y0); hold on;
AR2.EdgeColor   = 'k';
AR2.EdgeAlpha   = 1;
AR2.FaceColor   = 'r';
AR2.FaceAlpha   = 1;
AR              = area(x,yc.*y0); hold on;
AR.EdgeColor    = 'k';
AR.EdgeAlpha    = 1;
AR.FaceColor    = 'k';
AR.FaceAlpha    = 1;

plot(x,y0,'k','Linewidth',1); hold on;
axis tight
set(gcf,'color','w');
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fs)
set(gca,'FontName',fn)
set(gca,'Linewidth',lw)
