% Illustration particle filter vs Kalman filter, Figure 5.14

% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig
clear all
close all
clc
%% ------------------------------------------------------------------------ 
% Simulation univariate AR model (5.60)
% -------------------------------------------------------------------------
phi          = 0.7;                                                         % parameter
sw           = 0.2;                                                         % process noise std
sv           = 0.1;                                                         % measurement noise std
T            = 10;                                                          % Number of samples
x0           = 0;                                                           % initial state
x            = x0;                                                          % ...
y            = zeros(T,1);                                                  % Allocation simulated output
for k = 1:T                                                                 % simulate
    x        = phi*x + sw*randn(1,1);                                       % ...
    y(k)     = x + sv*randn(1,1);                                           % ...
end                                                                         % ...
%% ------------------------------------------------------------------------ 
% Kalman filter
% -------------------------------------------------------------------------
sx           = (sw^2)/(1-phi^2);                                            % initial state variance
Q            = sw^2;                                                        % process noise variance
R            = sv^2;                                                        % measurement noise variance
xKF          = zeros(1,T);                                                  % Allocation posterior state mean
PKF          = zeros(1,T);                                                  % Allocation posterior state covariance
xKF(1)       = x0 + sqrt(sx)*randn(1,1);                                    % Random state init
PKF(1)       = sx;                                                          % ...
for k = 1:T
    xKF(k)   = phi*xKF(k);
    PKF(k)   = phi*PKF(k)*phi' + Q;
    S        = PKF(k) + R;
    xKF(k)   = xKF(k) + PKF(k)/S*(y(k)-xKF(k));
    PKF(k)   = PKF(k) - PKF(k)/S*PKF(k);
end
%% ------------------------------------------------------------------------  
% Bootstrap Particle filter, Algorithm A.5
% -------------------------------------------------------------------------
plotheight   = 15;
plotwidth    = 15;
subplotsx    = 2;
subplotsy    = 5;
leftedge     = 0.5;
rightedge    = 0.5;
topedge      = 0.4;
bottomedge   = 0.6;
spacex       = 1;
spacey       = 0.9;
sub_pos      = subplot_pos(plotwidth,plotheight,leftedge,rightedge,...
                    bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
posY         = [5 5 4 4 3 3 2 2 1 1];                                       % trick to y position sub_pos
posX         = [1 2 1 2 1 2 1 2 1 2];                                       % trick to X position sub_pos
lw           = 1.5;
plotheight   = 15;
plotwidth    = 15;
fs           = 14;                                                          % fontSize
fw           = 'demi';                                                      % fontWeight
fn           = 'Times New Roman';                                           % fontName
figure('visible','on','units','normalized','outerposition',...
                                                    [0.5 0.04 0.5 0.95]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
dcolors      = linspecer(4,'qualitative');

N            = 1e6;                                                         % Number of particles
x            = xKF(1) + sqrt(sx)*randn(N,1);                                % 1. same initial conditions as the KF
for k = 1:T
    x        = phi.*x + sw.*randn(N,1);                                     % 2. p(x[k] | x[k-1])
    e        = repmat(y(k),N,1) - x;                                        % 3. Normalized importance weights
    wlog     = -0.5.*(log(2*pi) + log(R) + e./R.*e);                        % ...
    w        = exp(wlog - max(wlog));                                       % ...
    axes('position',sub_pos{posX(k),posY(k)});
    [fKF,xiKF] = ksdensity(xKF(k)+sqrt(PKF(k))*randn(1e6,1));
    [fPF,xiPF] = ksdensity(x);
    stem(x,w,'color',[211 211 211]./255,'marker','none','Linewidth',lw); hold on;
    plot(xiKF,fKF,'color',[255 0 0]./255,'Linewidth',lw); hold on; plot(xiPF,fPF,'color',dcolors(1,:),'Linewidth',lw); hold on;
    plot(y(k),max(fKF),'ko','Linewidth',1,'MarkerFaceColor','g'); hold on;
    text(0.1,0.9,['k = ',num2str(k)],'FontSize',fs,'FontWeight',fw,'FontName',fn,'Units','Normalized')    
    w        = w./sum(w);                                                   % ...
    idx      = sysresample(w);                                              % 4. Resample
    x        = x(idx);                                                      % ...
    [fPF,xiPF] = ksdensity(x);
    plot(xiPF,fPF,'color',[0 0 255]./255,'Linewidth',lw);
    axis tight
    set(gcf,'color','w');
    set(gca,'Box','off')
    set(gca,'FontWeight','demi')
    set(gca,'FontSize',fs)
    set(gca,'FontName',fn)
    set(gca,'Linewidth',lw)
end