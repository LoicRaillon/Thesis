% Plot figure 4.6, the figure adjustement depends on the screen
% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig

clear all
close all
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');
% -------------------------------------------------------------------------
% Import data for model M3
% -------------------------------------------------------------------------
load('M3_5alpha_chains.mat')
impM3                       = load('M3_5alpha_chains.mat');
stop                        = sys.MC;
start                       = 501;
chain1                      = reshape(p_ALL(1,:,:),[Nukwn,stop]);
chain2                      = reshape(p_ALL(2,:,:),[Nukwn,stop]);
chain3                      = reshape(p_ALL(3,:,:),[Nukwn,stop]);
chain5                      = reshape(p_ALL(5,:,:),[Nukwn,stop]);
chain6                      = reshape(p_ALL(6,:,:),[Nukwn,stop]);

% Inverse transform to constrained space
theta1_M3                   = zeros(Nukwn,stop);
theta2_M3                   = zeros(Nukwn,stop);
theta3_M3                   = zeros(Nukwn,stop);
theta4_M3                   = zeros(Nukwn,stop);
theta5_M3                   = zeros(Nukwn,stop);

for i = 1:stop
    theta1_M3(idx1,i)       = inv_transform(chain1(idx1,i),'LowUp',lb(idx1),ub(idx1));
    theta1_M3(idx2,i)       = inv_transform(chain1(idx2,i),'Log');
    theta1_M3(idx3,i)       = inv_transform(chain1(idx3,i),'LowUp',lb(idx3),ub(idx3));
    % --- %
    theta2_M3(idx1,i)       = inv_transform(chain2(idx1,i),'LowUp',lb(idx1),ub(idx1));
    theta2_M3(idx2,i)       = inv_transform(chain2(idx2,i),'Log');
    theta2_M3(idx3,i)       = inv_transform(chain2(idx3,i),'LowUp',lb(idx3),ub(idx3));
    % --- %
    theta3_M3(idx1,i)       = inv_transform(chain3(idx1,i),'LowUp',lb(idx1),ub(idx1));
    theta3_M3(idx2,i)       = inv_transform(chain3(idx2,i),'Log');
    theta3_M3(idx3,i)       = inv_transform(chain3(idx3,i),'LowUp',lb(idx3),ub(idx3));
    % --- %
    theta4_M3(idx1,i)       = inv_transform(chain5(idx1,i),'LowUp',lb(idx1),ub(idx1));
    theta4_M3(idx2,i)       = inv_transform(chain5(idx2,i),'Log');
    theta4_M3(idx3,i)       = inv_transform(chain5(idx3,i),'LowUp',lb(idx3),ub(idx3));
    % --- %
    theta5_M3(idx1,i)       = inv_transform(chain6(idx1,i),'LowUp',lb(idx1),ub(idx1));
    theta5_M3(idx2,i)       = inv_transform(chain6(idx2,i),'Log');
    theta5_M3(idx3,i)       = inv_transform(chain6(idx3,i),'LowUp',lb(idx3),ub(idx3));
end

sys.fun                     = @M3_5alpha;
thetaAll_M3                 = [theta1_M3(:,start:stop) theta2_M3(:,start:stop) theta3_M3(:,start:stop) theta4_M3(:,start:stop) theta5_M3(:,start:stop)];
Utheta                      = unique(thetaAll_M3','rows','stable')';
nUni                        = max(size(Utheta));
para                        = zeros(Np+2*Nx+Ny,1);
para(~fix,:)                = y(1);
para(fix,:)                 = Utheta(:,5899);
resM3                       = Residuals(para,sys);

% -------------------------------------------------------------------------
% Import data for model M4
% -------------------------------------------------------------------------
load('M4_5alpha_chain12.mat');
impM4_1                     = load('M4_5alpha_chain12.mat');
impM4_2                     = load('M4_5alpha_chain345.mat');
c1                          = sum(fix(1:Np));
c2                          = sum(fix(Np+1:Np+Nx+Ny));
c3                          = sum(fix(Np+Nx+Ny+1:Np+Nx+Ny+Nx));
idx11                       = 1:c1;
idx22                       = c1+1:c1+c2;
idx33                       = c1+c2+1:c1+c2+c3;

stop                        = sys.MC;
chain1                      = reshape(impM4_1.p_ALL(1,:,:),[Nukwn,stop]);
chain2                      = reshape(impM4_1.p_ALL(2,:,:),[Nukwn,stop]);
chain3                      = reshape(impM4_2.p_ALL(1,:,:),[Nukwn,stop]);
chain4                      = reshape(impM4_2.p_ALL(2,:,:),[Nukwn,stop]);
chain5                      = reshape(impM4_2.p_ALL(3,:,:),[Nukwn,stop]);

% Inverse transform to constrained space
theta1_M4                   = zeros(Nukwn,stop);
theta2_M4                   = zeros(Nukwn,stop);
theta3_M4                   = zeros(Nukwn,stop);
theta4_M4                   = zeros(Nukwn,stop);
theta5_M4                   = zeros(Nukwn,stop);

for i = 1:stop
    theta1_M4(idx11,i)      = inv_transform(chain1(idx11,i),'LowUp',lb(idx1),ub(idx1));
    theta1_M4(idx22,i)      = inv_transform(chain1(idx22,i),'Log');
    theta1_M4(idx33,i)      = inv_transform(chain1(idx33,i),'LowUp',lb(idx3),ub(idx3));
    % --- %
    theta2_M4(idx11,i)      = inv_transform(chain2(idx11,i),'LowUp',lb(idx1),ub(idx1));
    theta2_M4(idx22,i)      = inv_transform(chain2(idx22,i),'Log');
    theta2_M4(idx33,i)      = inv_transform(chain2(idx33,i),'LowUp',lb(idx3),ub(idx3));
    % --- %
    theta3_M4(idx11,i)      = inv_transform(chain3(idx11,i),'LowUp',lb(idx1),ub(idx1));
    theta3_M4(idx22,i)      = inv_transform(chain3(idx22,i),'Log');
    theta3_M4(idx33,i)      = inv_transform(chain3(idx33,i),'LowUp',lb(idx3),ub(idx3));
    % --- %
    theta4_M4(idx11,i)      = inv_transform(chain4(idx11,i),'LowUp',lb(idx1),ub(idx1));
    theta4_M4(idx22,i)      = inv_transform(chain4(idx22,i),'Log');
    theta4_M4(idx33,i)      = inv_transform(chain4(idx33,i),'LowUp',lb(idx3),ub(idx3));
    % --- %
    theta5_M4(idx11,i)      = inv_transform(chain5(idx11,i),'LowUp',lb(idx1),ub(idx1));
    theta5_M4(idx22,i)      = inv_transform(chain5(idx22,i),'Log');
    theta5_M4(idx33,i)      = inv_transform(chain5(idx33,i),'LowUp',lb(idx3),ub(idx3));
end

sys.fun                     = @M4_5alpha;
thetaAll_M4                 = [theta1_M4(:,start:stop) theta2_M4(:,start:stop) theta3_M4(:,start:stop) theta4_M4(:,start:stop) theta5_M4(:,start:stop)];
Utheta                      = unique(thetaAll_M4','rows','stable')';
nUni                        = max(size(Utheta));
para                        = zeros(Np+2*Nx+Ny,1);
para(~fix,:)                = [2.52e-3 ; 2.52e-4 ; 1e-6 ; y(1)];
para(fix,:)                 = Utheta(:,925);
resM4                       = Residuals(para,sys);

% -------------------------------------------------------------------------
% Plot configuration
% -------------------------------------------------------------------------
dcolors                     = linspecer(9,'qualitative');
% 1: red
% 2: blue
% 3: green
% 4: orange
% 5: yellow
% 6: brown
% 7: pink
% 8: grey
% 9: purple
lw                          = 2;
plotheight                  = 15;
plotwidth                   = 15;
fs                          = 14; % fontSize
fw                          = 'demi'; % fontWeight
fn                          = 'Times New Roman'; % fontName

% -------------------------------------------------------------------------
% Residuals analysis plot
% -------------------------------------------------------------------------
Nlags                       = 24*6;
N                           = length(y);
th                          = 1.96/sqrt(N);
ms                          = 3; % marker size
subplotsx                   = 2;
subplotsy                   = 4;
leftedge                    = 0.8;
rightedge                   = 0.4;
topedge                     = 0.3;
bottomedge                  = 0.95;
spacex                      = 1;
spacey                      = 0.9;
sub_pos                     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

figure('visible','on','units','normalized','outerposition',[0 0.0370 0.5000 0.9630]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------------------- %
axes('position',sub_pos{1,4});
% Autocorrelation %
s1                          = resM3;
s2                          = s1;
s1m                         = s1 - mean(s1);
s2m                         = s2 - mean(s2);
[cM3,~]                     = xcorr(s1m, s2m, Nlags, 'coeff');
s1                          = resM4;
s2                          = s1;
s1m                         = s1 - mean(s1);
s2m                         = s2 - mean(s2);
[cM4,lags]                  = xcorr(s1m, s2m, Nlags, 'coeff');
% remove lag at 0
h                           = stem(lags(Nlags+1:end),[cM3(Nlags+1:end) cM4(Nlags+1:end)],'filled','linewidth',lw-1);
h(1).Color                  = dcolors(2,:);
h(1).MarkerSize             = ms;
h(2).Color                  = dcolors(6,:);
h(2).MarkerSize             = ms;
line([-2 Nlags+2],[th th],'color',[1 0 0],'linewidth',1)
line([-2 Nlags+2],[-th -th],'color',[1 0 0],'linewidth',1)
axis tight
xlim([-2 Nlags+2])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xlabel('lags','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.3,0.75,'\textbf{ACF($\bar{e}$)}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------------------- %
axes('position',sub_pos{2,4});
cpgram(resM3,1,dcolors(2,:)); hold on;
cpgram(resM4,1,dcolors(6,:)); hold on;
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
text(0.3,0.75,'\textbf{CP($\bar{e}$)}','FontSize',fs,'FontWeight','demi','FontName','Cambria','Units','Normalized','interpreter','latex')
% --------------------------------------- %
axes('position',sub_pos{1,3});
s2                          = u(1,:);
s2m                         = s2 - mean(s2);
s1                          = resM3;
s1m                         = s1 - mean(s1);
[cM3,~]                     = xcorr(s1m, s2m, Nlags, 'coeff');
s1                          = resM4;
s1m                         = s1 - mean(s1);
[cM4,lags]                  = xcorr(s1m, s2m, Nlags, 'coeff');
% remove lag at 0
h                           = stem(lags(Nlags+1:end),[cM3(Nlags+1:end) cM4(Nlags+1:end)],'filled','linewidth',lw-1);
h(1).Color                  = dcolors(2,:);
h(1).MarkerSize             = ms;
h(2).Color                  = dcolors(6,:);
h(2).MarkerSize             = ms;
line([-2 Nlags+2],[th th],'color',[1 0 0],'linewidth',1)
line([-2 Nlags+2],[-th -th],'color',[1 0 0],'linewidth',1)
axis tight
xlim([-2 Nlags+2])
ylim([-0.07 0.09])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xlabel('lags','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.3,0.85,'\textbf{CCF($\bar{e},T_o$)}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------------------- %
axes('position',sub_pos{2,3});
s2                          = u(2,:);
s2m                         = s2 - mean(s2);
s1                          = resM3;
s1m                         = s1 - mean(s1);
[cM3,~]                     = xcorr(s1m, s2m, Nlags, 'coeff');
s1                          = resM4;
s1m                         = s1 - mean(s1);
[cM4,lags]                  = xcorr(s1m, s2m, Nlags, 'coeff');
% remove lag at 0
h                           = stem(lags(Nlags+1:end),[cM3(Nlags+1:end) cM4(Nlags+1:end)],'filled','linewidth',lw-1);
h(1).Color                  = dcolors(2,:);
h(1).MarkerSize             = ms;
h(2).Color                  = dcolors(6,:);
h(2).MarkerSize             = ms;
line([-2 Nlags+2],[th th],'color',[1 0 0],'linewidth',1)
line([-2 Nlags+2],[-th -th],'color',[1 0 0],'linewidth',1)
axis tight
xlim([-2 Nlags+2])
ylim([-0.07 0.09])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xlabel('lags','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.3,0.85,'\textbf{CCF($\bar{e},T_n$)}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------------------- %
axes('position',sub_pos{1,2});
s2                          = u(3,:);
s2m                         = s2 - mean(s2);
s1                          = resM3;
s1m                         = s1 - mean(s1);
[cM3,~]                     = xcorr(s1m, s2m, Nlags, 'coeff');
s1                          = resM4;
s1m                         = s1 - mean(s1);
[cM4,lags]                  = xcorr(s1m, s2m, Nlags, 'coeff');
% remove lag at 0
h                           = stem(lags(Nlags+1:end),[cM3(Nlags+1:end) cM4(Nlags+1:end)],'filled','linewidth',lw-1);
h(1).Color                  = dcolors(2,:);
h(1).MarkerSize             = ms;
h(2).Color                  = dcolors(6,:);
h(2).MarkerSize             = ms;
line([-2 Nlags+2],[th th],'color',[1 0 0],'linewidth',1)
line([-2 Nlags+2],[-th -th],'color',[1 0 0],'linewidth',1)
axis tight
xlim([-2 Nlags+2])
ylim([-0.07 0.09])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xlabel('lags','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.3,0.85,'\textbf{CCF($\bar{e},\dot{Q}_{W}$)}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------------------- %
axes('position',sub_pos{2,2});
s2                          = u(4,:);
s2m                         = s2 - mean(s2);
s1                          = resM3;
s1m                         = s1 - mean(s1);
[cM3,~]                     = xcorr(s1m, s2m, Nlags, 'coeff');
s1                          = resM4;
s1m                         = s1 - mean(s1);
[cM4,lags]                  = xcorr(s1m, s2m, Nlags, 'coeff');
% remove lag at 0
h                           = stem(lags(Nlags+1:end),[cM3(Nlags+1:end) cM4(Nlags+1:end)],'filled','linewidth',lw-1);
h(1).Color                  = dcolors(2,:);
h(1).MarkerSize             = ms;
h(2).Color                  = dcolors(6,:);
h(2).MarkerSize             = ms;
line([-2 Nlags+2],[th th],'color',[1 0 0],'linewidth',1)
line([-2 Nlags+2],[-th -th],'color',[1 0 0],'linewidth',1)
axis tight
xlim([-2 Nlags+2])
ylim([-0.07 0.09])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xlabel('lags','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.3,0.85,'\textbf{CCF($\bar{e},\dot{Q}_{E}$)}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------------------- %
axes('position',sub_pos{1,1});
s2                          = u(5,:);
s2m                         = s2 - mean(s2);
s1                          = resM3;
s1m                         = s1 - mean(s1);
[cM3,~]                     = xcorr(s1m, s2m, Nlags, 'coeff');
s1                          = resM4;
s1m                         = s1 - mean(s1);
[cM4,lags]                  = xcorr(s1m, s2m, Nlags, 'coeff');
% remove lag at 0
h                           = stem(lags(Nlags+1:end),[cM3(Nlags+1:end) cM4(Nlags+1:end)],'filled','linewidth',lw-1);
h(1).Color                  = dcolors(2,:);
h(1).MarkerSize             = ms;
h(2).Color                  = dcolors(6,:);
h(2).MarkerSize             = ms;
line([-2 Nlags+2],[th th],'color',[1 0 0],'linewidth',1)
line([-2 Nlags+2],[-th -th],'color',[1 0 0],'linewidth',1)
axis tight
xlim([-2 Nlags+2])
ylim([-0.07 0.09])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xlabel('lags','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.3,0.85,'\textbf{CCF($\bar{e},\dot{Q}_{S}$)}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------------------- %
axes('position',sub_pos{2,1});
s2                          = u(6,:);
s2m                         = s2 - mean(s2);
s1                          = resM3;
s1m                         = s1 - mean(s1);
[cM3,~]                     = xcorr(s1m, s2m, Nlags, 'coeff');
s1                          = resM4;
s1m                         = s1 - mean(s1);
[cM4,lags]                  = xcorr(s1m, s2m, Nlags, 'coeff');
% remove lag at 0
h                           = stem(lags(Nlags+1:end),[cM3(Nlags+1:end) cM4(Nlags+1:end)],'filled','linewidth',lw-1);
h(1).Color                  = dcolors(2,:);
h(1).MarkerSize             = ms;
h(2).Color                  = dcolors(6,:);
h(2).MarkerSize             = ms;
line([-2 Nlags+2],[th th],'color',[1 0 0],'linewidth',1)
line([-2 Nlags+2],[-th -th],'color',[1 0 0],'linewidth',1)
axis tight
xlim([-2 Nlags+2])
ylim([-0.07 0.09])
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
xlabel('lags','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.3,0.85,'\textbf{CCF($\bar{e},\dot{Q}_h$)}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
