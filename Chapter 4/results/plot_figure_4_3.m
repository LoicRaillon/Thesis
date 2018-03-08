% Plot figure 4.3, the figure adjustement depends on the screen
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
para                        = zeros(Np+2*Nx+Ny,1);
para(~fix,:)                = [2.52e-3 ; 2.52e-4 ; 1e-6 ; y(1)];
para(fix,:)                 = Utheta(:,925);
resM4                       = Residuals(para,sys);

% -------------------------------------------------------------------------
% Plot configuration
% -------------------------------------------------------------------------
dcolors                     = linspecer(9,'qualitative');
% 0: black --> Ts
% 1: red
% 2: blue --> M3
% 3: green --> Qe
% 4: orange --> Qw
% 5: yellow
% 6: brown --> M4
% 7: pink 
% 8: grey --> Qh
% 9: purple --> Qs
lw                          = 2;
plotheight                  = 15;
plotwidth                   = 15;
fs                          = 14; % fontSize
fw                          = 'demi'; % fontWeight
fn                          = 'Times New Roman'; % fontName

cutp                        = 2890:6441;
Tp                          = length(cutp);
y                           = Tzone1_all(cutp);
u                           = [Text(cutp) Tzone2_all(cutp) Sr_west(cutp) Sr_east(cutp) Sr_south(cutp) Qhvac1(cutp)]';
startDate                   = datenum('29-04-2014 01:50:00', 'dd-mm-yy HH:MM:SS');
endDate                     = datenum('23-05-2014 17:30:00', 'dd-mm-yy HH:MM:SS'); 
xData                       = linspace(startDate,endDate,Tp);

subplotsx                   = 1;
subplotsy                   = 6;
leftedge                    = 1;
rightedge                   = 0.6;
topedge                     = 0.3;
bottomedge                  = 0.7;
spacex                      = 1.3;
spacey                      = 0.9;
sub_pos                     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

figure('visible','on','units','normalized','outerposition',[0 0.0370 0.5000 0.9630]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------------------- %
axes('position',sub_pos{1,6});
plot(xData,y,'k','linewidth',lw-1); hold on;
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
bis = get(gca);
plot([xData(2111) xData(2111)],bis.YLim,'--k','linewidth',lw-1);
ylabel('°C','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.025,0.9,'\textbf{$T_s$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex');
% --------------------------------------- %
axes('position',sub_pos{1,5});
ga = [resM3;zeros(1441,1)];
plot(xData,ga,'Color',dcolors(2,:),'linewidth',lw-1); hold on;
plot(xData(2112:3552),ga(2112:3552),'w','linewidth',lw-1); hold on;
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
bis = get(gca);
plot([xData(2111) xData(2111)],bis.YLim,'--k','linewidth',lw-1);
ylabel('°C','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.025,0.9,'\textbf{$\bar{e}$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex','Color','k');
% --------------------------------------- %
axes('position',sub_pos{1,4});
plot(xData,[resM4;zeros(1441,1)],'Color',dcolors(6,:),'linewidth',lw-1); hold on;
plot(xData(2112:3552),ga(2112:3552),'w','linewidth',lw-1); hold on;
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
bis = get(gca);
plot([xData(2111) xData(2111)],bis.YLim,'--k','linewidth',lw-1);
ylabel('°C','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.025,0.9,'\textbf{$\bar{e}$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex','Color','k');
% --------------------------------------- %
axes('position',sub_pos{1,3});
plot(xData,u(1,:),'Color','k','linewidth',lw-1); hold on; %dcolors(1,:)
plot(xData,u(2,:),'Color','k','linewidth',lw-1); hold on; % dcolors(7,:)
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
bis = get(gca);
plot([xData(2111) xData(2111)],bis.YLim,'--k','linewidth',lw-1);
ylabel('°C','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.025,0.6,'\textbf{$T_o$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex','Color','k');
text(0.025,0.95,'\textbf{$T_n$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex','Color','k');
% --------------------------------------- %
axes('position',sub_pos{1,2});
plot(xData,u(3,:),'Color',dcolors(4,:),'linewidth',lw-1); hold on; 
plot(xData,u(4,:),'Color',dcolors(3,:),'linewidth',lw-1); hold on;
plot(xData,u(5,:),'Color',dcolors(9,:),'linewidth',lw-1); hold on;
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
x1 = get(gca,'Ytick');
set(gca,'YTickLabel',x1./1000)
bis = get(gca);
plot([xData(2111) xData(2111)],bis.YLim,'--k','linewidth',lw-1);
ylabel('kW','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.025,1.05,'\textbf{$\dot{Q}_{W}$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex','Color',dcolors(4,:));
text(0.065,1.05,'\textbf{$\dot{Q}_{E}$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex','Color',dcolors(3,:));
text(0.105,1.05,'\textbf{$\dot{Q}_{S}$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex','Color',dcolors(9,:));
% --------------------------------------- %
axes('position',sub_pos{1,1});
plot(xData,u(6,:),'Color','k','linewidth',lw-1); hold on; %dcolors(8,:)
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
x1 = get(gca,'Ytick');
set(gca,'YTickLabel',x1./1000)
bis = get(gca);
plot([xData(2111) xData(2111)],bis.YLim,'--k','linewidth',lw-1);
ylabel('kW','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(0.025,0.95,'\textbf{$\dot{Q}_h$}','FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex','Color','k');
