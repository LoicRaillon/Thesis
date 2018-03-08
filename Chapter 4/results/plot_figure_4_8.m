% Figure 4.8: Illustration of the penalty function 
% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - export_fig: Export figure: available at: https://github.com/altmany/export_fig

clear all
close all
clc
addpath(userpath, 'utilities');
%% ------------------------------------------------------------------------
% Setup signal
% -------------------------------------------------------------------------
N               = 1e6;                              % Number of samples
Nl              = 3;                                % Number of lambda (penalty parameter)
dcolors         = linspecer(Nl+1,'qualitative');    % Distinguishable color
lambda          = linspace(1e-6,1e-1,Nl);           % Nl values between 1e-6 and 1e-1
ub              = 0.8;                              % Upper bound
lb              = -0.5;                             % Lower bound
x               = linspace(lb,ub,N);                % 1e6 values between lb and ub
f               = x.^2;                             % Square function 
fpen            = zeros(Nl,N);                      % Allocation penalysed square function 

%% ------------------------------------------------------------------------
% Setup figure
% -------------------------------------------------------------------------
plotheight      = 15;
plotwidth       = 15;
lw              = 1.5;                              % LineWidth
fs              = 14;                               % fontSize
fw              = 'demi';                           % fontWeight
fn              = 'Times New Roman';                % fontName
figure('visible','on','units','normalized','outerposition',[0 0.5 0.4 0.5]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');

%% ------------------------------------------------------------------------
% PLot f(x) = x² and fpen(x) = x² + pen(x)
% -------------------------------------------------------------------------
for i = 1:Nl
   fpen(i,:)    = x.^2 + lambda(i).*( abs(lb)./(x-lb) + abs(ub)./(ub-x) ); 
   bis          = x(fpen(i,:)==min(fpen(i,:)));
   plot(x,fpen(i,:),'Color',dcolors(i+1,:),'linewidth',lw); hold on;
   plot([bis bis],[0 1],'Color',dcolors(i+1,:),'linewidth',lw); hold on
end
plot(x,f,'Color',dcolors(1,:),'linewidth',lw,'LineStyle','--'); hold on;
xlim([lb ub])
ylim([0 1])
set(gcf,'color','w');
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fs)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',lw)
