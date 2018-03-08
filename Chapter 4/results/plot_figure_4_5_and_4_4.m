% Plot figure 4.5, the figure adjustement depends on the screen
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
theta1_M4([6 7],:)          = theta1_M4([7 6],:);
theta2_M4([6 7],:)          = theta2_M4([7 6],:);
theta3_M4([6 7],:)          = theta3_M4([7 6],:);
theta4_M4([6 7],:)          = theta4_M4([7 6],:);
theta5_M4([6 7],:)          = theta5_M4([7 6],:);

% -------------------------------------------------------------------------
% Plot configuration
% -------------------------------------------------------------------------
dcolors             = linspecer(9,'qualitative');
% 1: red
% 2: blue
% 3: green
% 4: orange 
% 5: yellow
% 6: brown
% 7: pink
% 8: grey
% 9: purple
lw                  = 2;
plotheight          = 15;
plotwidth           = 15;
fs                  = 14; % fontSize
fw                  = 'demi'; % fontWeight
fn                  = 'Times New Roman'; % fontName
% Cm and Ci inverted in M3 and M4
VarNames            = {'\textbf{$R_o$}','\textbf{$R_i$}','\textbf{$R_m$}','\textbf{$R_z$}', ...
    '\textbf{$C_w$}','\textbf{$C_m$}','\textbf{$C_i$}', ...
    '\textbf{$\alpha_{w_{W}}$}','\textbf{$\alpha_{w_{S}}$}', ...
    '\textbf{$\alpha_{i_{W}}$}','\textbf{$\alpha_{i_{E}}$}','\textbf{$\alpha_{i_{S}}$}', ...
    '\textbf{$\sigma_{w_{11}}$}','\textbf{$\sigma_{w_{22}}$}','\textbf{$\sigma_{w_{33}}$}', ...
    '\textbf{$\sigma_v$}','\textbf{$x_{w_{0}}$}','\textbf{$x_{m_{0}}$}','\textbf{$x_{i_{0}}$}'};
posY                = [7 6 5 4 3 2 1]; % trick to y position sub_pos

% -------------------------------------------------------------------------
% Posterior distributions plot
% -------------------------------------------------------------------------
start               = 501;
thetaAll_M3         = [theta1_M3(:,start:stop) theta2_M3(:,start:stop) theta3_M3(:,start:stop) theta4_M3(:,start:stop) theta5_M3(:,start:stop)];
thetaAll_M4         = [theta1_M4(:,start:stop) theta2_M4(:,start:stop) theta3_M4(:,start:stop) theta4_M4(:,start:stop) theta5_M4(:,start:stop)];

subplotsx           = 3;
subplotsy           = 7;
leftedge            = 0.8;
rightedge           = 0.8;
topedge             = 0.5;
bottomedge          = 0.6;
spacex              = 1.3;
spacey              = 0.9;
alpha               = 0.25;
sub_pos             = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

figure('visible','on','units','normalized','outerposition',[0 0.0370 0.5000 0.9630]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');

k                   = 1;
for i = 1:subplotsy
    for j = 1:subplotsx
        if k == 20
            break
        end
        axes('position',sub_pos{j,posY(i)});
        histogram(theta1_M4(k,start:stop),'FaceColor',dcolors(1,:),'FaceAlpha',alpha,'EdgeColor',dcolors(1,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
        histogram(theta2_M4(k,start:stop),'FaceColor',dcolors(3,:),'FaceAlpha',alpha,'EdgeColor',dcolors(2,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
        histogram(theta3_M4(k,start:stop),'FaceColor',dcolors(4,:),'FaceAlpha',alpha,'EdgeColor',dcolors(3,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
        histogram(theta4_M4(k,start:stop),'FaceColor',dcolors(7,:),'FaceAlpha',alpha,'EdgeColor',dcolors(4,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
        histogram(theta5_M4(k,start:stop),'FaceColor',dcolors(8,:),'FaceAlpha',alpha,'EdgeColor',dcolors(5,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
        [f,xi]  = ksdensity(thetaAll_M4(k,:));
        expectM4(k) = xi(f == max(f));
        plot(xi,f,'Color',dcolors(6,:),'Linewidth',lw)
        if k < 19
            histogram(theta1_M3(k,start:stop),'FaceColor',dcolors(1,:),'FaceAlpha',alpha,'EdgeColor',dcolors(1,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
            histogram(theta2_M3(k,start:stop),'FaceColor',dcolors(3,:),'FaceAlpha',alpha,'EdgeColor',dcolors(2,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
            histogram(theta3_M3(k,start:stop),'FaceColor',dcolors(4,:),'FaceAlpha',alpha,'EdgeColor',dcolors(3,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
            histogram(theta4_M3(k,start:stop),'FaceColor',dcolors(7,:),'FaceAlpha',alpha,'EdgeColor',dcolors(4,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
            histogram(theta5_M3(k,start:stop),'FaceColor',dcolors(8,:),'FaceAlpha',alpha,'EdgeColor',dcolors(5,:),'EdgeAlpha',0,'Normalization','pdf'); hold on
            [f,xi]  = ksdensity(thetaAll_M3(k,:));
            expectM3(k) = xi(f == max(f));
            plot(xi,f,'Color',dcolors(2,:),'Linewidth',lw)
        end
        axis tight
        set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
        text(0.95,0.9,VarNames{k},'FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
        xx          = get(gca,'Xtick');
        yy          = get(gca,'Ytick');
        nx          = max(floor(log(abs(xx))./log(10)));
        ny          = max(floor(log(abs(yy))./log(10)));
        if nx <= -1
            set(gca,'XTickLabel',xx.*(10^-nx))
            text(1.02,0.05,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
        end
        if nx >= 1
            set(gca,'XTickLabel',xx./(10^nx))
            text(1.02,0.05,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
        end
        if ny >= 1
            set(gca,'YTickLabel',yy./(10^ny))
            text(-0.04,1.18,['10^{',num2str(ny),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
        end
        k           = k + 1;
    end
end

return
% -------------------------------------------------------------------------
% Trace plot M3
% -------------------------------------------------------------------------
posY                = [6 5 4 3 2 1]; % trick to y position sub_pos
subplotsx           = 3;
subplotsy           = 6;
leftedge            = 0.8;
rightedge           = 0.8;
topedge             = 0.5;
bottomedge          = 0.6;
spacex              = 1.3;
spacey              = 0.9;
alpha               = 0.25;
sub_pos             = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

figure('visible','on','units','normalized','outerposition',[0 0.0370 0.5000 0.9630]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');

k                   = 1;
for i = 1:subplotsy
    for j = 1:subplotsx
        axes('position',sub_pos{j,posY(i)});
        plot(theta1_M3(k,1:2e2),'Color',dcolors(1,:),'Linewidth',lw); hold on
        plot(theta2_M3(k,1:2e2),'Color',dcolors(3,:),'Linewidth',lw); hold on
        plot(theta3_M3(k,1:2e2),'Color',dcolors(4,:),'Linewidth',lw); hold on
        plot(theta4_M3(k,1:2e2),'Color',dcolors(7,:),'Linewidth',lw); hold on
        plot(theta5_M3(k,1:2e2),'Color',dcolors(8,:),'Linewidth',lw); hold on
        axis tight
        set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
        text(0.95,1,VarNames{k},'FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
        xx          = get(gca,'Xtick');
        yy          = get(gca,'Ytick');
        nx          = max(floor(log(abs(xx))./log(10)));
        ny          = max(floor(log(abs(yy))./log(10)));
        if nx <= -1
            set(gca,'XTickLabel',xx.*(10^-nx))
            text(1.02,0.05,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
        end
        if nx >= 1
            set(gca,'XTickLabel',xx./(10^nx))
            text(1.02,0.05,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
        end
        if ny <= -1
            set(gca,'YTickLabel',yy.*(10^-ny))
            text(-0.04,1.18,['10^{',num2str(ny),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
        end
        if ny >= 1
            set(gca,'YTickLabel',yy./(10^ny))
            text(-0.04,1.18,['10^{',num2str(ny),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
        end
        k           = k + 1;
    end
end
