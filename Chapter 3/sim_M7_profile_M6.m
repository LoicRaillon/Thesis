% Simulation model detailed M7 and profile model M6
% Plot
%   - Figure 3.12: Inputs data
%   - Figure 3.13: Output and frequency spectrum
%   - Figure 3.14: Profile likelihood model M6
% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig

clear all;
close all;
clc
addpath(pwd, 'models');
addpath(pwd, 'utilities');
%% ------------------------------------------------------------------------
% Import data | Inputs and Outputs definition
% Import Twin Houses experiment data available at 
% https://pure.strath.ac.uk/portal/en/datasets/twin-houses-empirical-validation-dataset-experiment-2(94559779-e781-4318-8842-80a2b1201668).html
% The data of the Twin Houses experiment are used for the simulated examples
% -------------------------------------------------------------------------
dt                          = 10*60;                                        % sampling time [s]
imp                         = importdata('Twin_house_exp2_10min.xlsx');     % import data
cut                         = 505:1368;                                     % part of the data
T                           = length(cut);                                  % length of selected part
Text                        = imp.data(:,54);                               % Outdoor air temperature [°C]
Sr_north                    = imp.data(:,58);                               % solar radiation measured on the north vertical wall [W/m²]
Sr_east                     = imp.data(:,59);                               % solar radiation measured on the east vertical wall [W/m²]
Sr_south                    = imp.data(:,60);                               % solar radiation measured on the south vertical wall [W/m²]
Sr_west                     = imp.data(:,61);                               % solar radiation measured on the west vertical wall [W/m²]
%% ------------------------------------------------------------------------
% Simulate detailed model M7
% -------------------------------------------------------------------------
width                       = [5 7.5 7.5 5];                                % width of walls [North East West South] [m]
h                           = 2;                                            % heigth of walls [m]
Sw                          = width*h;                                      % Surfaces [m²]
S                           = sum(Sw);                                      % total surface [m²]
e                           = 0.2;                                          % Wall thickness [m]
hco                         = 25;                                           % outdoor convective coefficient [W/(m².K)]
hci                         = 8;                                            % indoor convective coefficient [W/(m².K)]
lambda                      = 2;                                            % Thermal conductivity [W/(m.K)]
rho                         = 2400;                                         % density [kg.m^-3]
c                           = 1000;                                         % specific heat capacity [J/(Kg.K)]
% Rso: Outdoor convective resistance [K/W]
% Rw: Thermal resistance of the concrete wall [K/W] 
% Rwo = Rwi = Rw/2: symmetric walls
% Rsi: Indoor convective resistance [K/W]
% Cw: Thermal capacity of the concrete wall [K/W]
% --------------------------North---------------------------------------- %
RsoN                        = 1/(hco*Sw(1));
RwN                         = e/(Sw(1)*lambda);
RwoN                        = RwN/2;
RwiN                        = RwN/2;
RsiN                        = 1/(hci*Sw(1));
CwN                         = rho*Sw(1)*e*c;
% --------------------------East----------------------------------------- %
RsoE                        = 1/(hco*Sw(2));
RwE                         = e/(Sw(2)*lambda);
RwoE                        = RwE/2;
RwiE                        = RwE/2;
RsiE                        = 1/(hci*Sw(2));
CwE                         = rho*Sw(2)*e*c;
% --------------------------West----------------------------------------- %
RsoW                        = 1/(hco*Sw(3));
RwW                         = e/(Sw(3)*lambda);
RwoW                        = RwW/2;
RwiW                        = RwW/2;
RsiW                        = 1/(hci*Sw(3));
CwW                         = rho*Sw(3)*e*c;
% --------------------------South---------------------------------------- %
RsoS                        = 1/(hco*Sw(4));
RwS                         = e/(Sw(4)*lambda);
RwoS                        = RwS/2;
RwiS                        = RwS/2;
RsiS                        = 1/(hci*Sw(4));
CwS                         = rho*Sw(4)*e*c;
Ci                          = 5.621e5;                                      % selected s.t the second time constant is 30 min 
% ----------------------------------------------------------------------- %
R                           = [RsoN RwoN RwiN RsiN RsoE RwoE RwiE RsiE ...  % Thermal resistances M7
                                RsoW RwoW RwiW RsiW RsoS RwoS RwiS RsiS];   % ...
C                           = [zeros(1,8) CwN CwE CwW CwS Ci];              % Thermal capacities M7
alpha                       = 0.6;                                          % Wall absorptivity
[A,B]                       = model_detailed(R,C);                          % Evaluate M7
C                           = [0 0 0 0 1];                                  % Output matrix
Nx                          = 5;                                            % Number of states
Ny                          = 1;                                            % Number of outputs
F                           = expm(A*dt);                                   % discretization
G                           = A\(F-eye(Nx))*B;                              % ...
Qstd                        = eye(Nx)*1e-2;                                 % std process noise
Rstd                        = eye(Ny)*1e-2;                                 % std measurement noise
xw0N                        = 15;                                           % Initial states
xw0E                        = 15;                                           % ...
xw0W                        = 15;                                           % ...
xw0S                        = 15;                                           % ...
xi0                         = 20;                                           % ...
x0                          = [xw0N xw0E xw0W xw0S xi0]';                   % ...
To                          = Text(cut);                                    % Inputs
Qo                          = alpha.*Sw.*[Sr_north(cut) Sr_east(cut) ...    % ...
                                            Sr_west(cut) Sr_south(cut)];    % ...                                     
%% ------------------------------------------------------------------------
% ROLBS signal: y1 --> Persistent excitation
% time constants from 20min to 1h --> fast and
% from 5h to 15h --> slow
% -------------------------------------------------------------------------
x                           = x0;                                           % Pass initial states
y1                          = zeros(Ny,T);                                  % Allocation simulated output
n                           = 3;
slow                        = ROLBS(n,666,30);
fast                        = ROLBS(n,666,2);
bis                         = [slow circshift(slow,round(n/3)) circshift(slow,round(2*n/3))];
bis2                        = [slow circshift(fast,round(n/3)) circshift(fast,round(2*n/3))];
Qh                          = 3500.*[bis bis2]';
Qh                          = Qh(1:864);                                    % 6*24*6
u1                          = [To Qo Qh]';                                  % input vector
for k = 1:T                                                                 % simulate with ROLBS signal
    y1(k)                   = C*x + Rstd*randn(Ny,1);                       % ...
    x                       = F*x + G*u1(:,k) + Qstd*randn(Nx,1);           % ...
end
%% ------------------------------------------------------------------------
% controled signal: y2 --> low excitation
% Kalman filter with unknown input
% for details see: Input and State Estimation for Linear Systems: A Least
% Squares Estimation Approach, Shuwen Pan, Hongye Su, Hong Wang, Senior
% Member, IEEE, Jian Chu and Renquan Lu
% Proceedings of the 7th Asian Control Conference, Hong Kong, China, August
% 27-29, 2009
% -------------------------------------------------------------------------
x                           = x0;
Q                           = Qstd.^2;
R                           = Rstd.^2;
P                           = eye(Nx)*0.01;
Gd                          = G(:,6);
Iy                          = eye(Ny);
Ix                          = eye(Nx);
dS                          = zeros(T,1);
d                           = 1.5e3;
dS(1)                       = d;
xS2                         = zeros(Nx,T);
xS2(:,1)                    = x;
y                           = ones(T,1)*20;
u                           = [To Qo zeros(T,1)]';
for k = 2:T
    xbar                    = F*x + G*u(:,k-1);
    x                       = xbar + Gd*d;
    P                       = F*P*F' + Q;
    Se                      = R + C*P*C';
    K                       = P*C'/Se;
    S                       = Gd'*C'/R*(Iy-C*K)*C*Gd;
    v1                      = y(k)-C*x;
    x                       = x + K*v1;
    v2                      = y(k)-C*xbar;
    d                       = S\Gd'*C'/R*(Iy-C*K)*v2;
    P                       = (Ix-K*C)*(P+Gd/S*Gd'*(Ix-K*C)');
    dS(k)                   = d;
    xS2(:,k)                = x;
end
Qh2                         = dS;
x                           = x0;
y2                          = zeros(Ny,T);
u2                          = [To Qo Qh2]';
for k = 1:T
    y2(k)                   = C*x + Rstd*randn(Ny,1);
    x                       = F*x + G*u2(:,k) + Qstd*randn(Nx,1);
end
%% ------------------------------------------------------------------------
% Set structure for profile loglikelihood model M6
% -------------------------------------------------------------------------
Np                          = 5;                                            % Number of R,C parameters
Nx                          = 2;                                            % Number of states
Ny                          = 1;                                            % Number of outputs
Nu                          = 3;                                            % Number of inputs
fix                         = [true(Np,1) ; false(2*Nx+Ny,1)];              % Logical indexing of fixed parameters
sys.dim                     = [Np Nx Ny Nu];                                % vector of dimensions
sys.fix                     = fix;                                          % Passing to the structure
sys.fun                     = @M6;                                          % Model M6
sys.Psqrt                   = diag([0.5 0.5]);                              % Initial std of states
sys.sc                      = 1e8;                                          % Scaling of thermal capacities
sys.dt                      = dt;                                           % Sampling time
sys.C                       = [0 1];                                        % output matrix

Rso                         = 1/(1/RsoN + 1/RsoE + 1/RsoW + 1/RsoS);        % Compute parameter values for model M6
Rsi                         = 1/(1/RsiN + 1/RsiE + 1/RsiW + 1/RsiS);        % ...
Rw                          = 1/(1/RwN + 1/RwE + 1/RwW + 1/RwS);            % ...
Rwo                         = Rw/2;                                         % ...
Rwi                         = Rw/2;                                         % ...
Ro                          = Rwo + Rso;                                    % ...
Ri                          = Rwi + Rsi;                                    % ...
Cw                          = CwN + CwE + CwW + CwS;                        % ...
aW                          = ((alpha*RsoN)/(RwoN+RsoN)+...                 % ...
                                    (alpha*RsoE)/(RwoE+RsoE))/2;            % ...
sigw1                       = 1e-2;                                         % std process noise
sigw2                       = 1e-2;                                         % ...
sigv                        = 1e-2;                                         % std measurement noise
xw0                         = 15;                                           % initial states
xi0                         = 20;                                           % ...
%                               lb     p     ub
gather(:,1)                 = [1e-4 ; Ro        ; 1];
gather(:,2)                 = [1e-4 ; Ri        ; 1];
gather(:,3)                 = [1e-4 ; Cw/sys.sc ; 10];
gather(:,4)                 = [1e-5 ; Ci/sys.sc ; 1];
gather(:,5)                 = [1e-4 ; aW        ; 10];
gather(:,6)                 = [1e-8 ; sigw1     ; 1];
gather(:,7)                 = [1e-8 ; sigw2     ; 1];
gather(:,8)                 = [1e-8 ; sigv      ; 1];
gather(:,9)                 = [10   ; xw0       ; 40];
gather(:,10)                = [10   ; xi0       ; 40];
p                           = gather(2,:)';                                 % Parameter vector
lb                          = gather(1,:)';                                 % lower bounds
ub                          = gather(3,:)';                                 % upper bounds
sys.lb                      = lb;                                           % Passing to the structure
sys.ub                      = ub;                                           % ...

Qo2                         = sum(Sw.*[Sr_north(cut) Sr_east(cut) ...       % sum of vertical solar radiation
                                        Sr_west(cut) Sr_south(cut)],2);     % ...
u3                          = [To Qo2 Qh]';                                 % input vector ROLBS
u4                          = [To Qo2 Qh2]';                                % input vector low excitation
options                     = optimoptions(@fminunc,...
                                    'Display','none',...
                                    'Algorithm','quasi-newton',...
                                    'FiniteDifferenceType','central',...
                                    'HessUpdate','bfgs',...
                                    'MaxIterations',1e3,...
                                    'MaxFunctionEvaluations',1e3,...
                                    'OptimalityTolerance',1e-9,...
                                    'UseParallel',true,...
                                    'StepTolerance',1e-9);
%% ------------------------------------------------------------------------
% Start profile loglikelihood model M6, low excitation 
% -------------------------------------------------------------------------
sys.u                       = u4;                                           % Find MLE
sys.y                       = y2;                                           % ...
x0                          = p(Np+Nx+Ny+1:end);                            % ...
eta                         = transform(p(fix),'LowUp',lb(fix),ub(fix));    % ...
fun                         = @(x) LogLikPen([x(1:5);p(6:8);x0],sys);       % ...
[eta_MLE,~,~,~,~,~]         = fminunc(fun,eta,options);                     % ...
theta_MLE2                  = inv_transform(eta_MLE,'LowUp',...             % ...
                                                lb(fix),ub(fix));           % ...
Npoints                     = 200;                                          % Number of grid points
prct                        = [0.6 0.3 0.3 0.4 0.6];                        % percentage of variation of the parameters
grid2                       = zeros(Np,Npoints+2);                          % Allocation grid
Ngrid                       = max(size(grid2));                             % size grid
PL2                         = zeros(Np,Ngrid);                              % Allocation profile

for j = 1:Np
    disp(['PL2: ',num2str(j),' / ',num2str(Np)])
    fix2                    = fix;                                          % Fix profile parameter
    fix2(j)                 = false;                                        % ...
    sys.fix                 = fix2;                                         % ...
    eta                     = transform(p(fix2),'LowUp',lb(fix2),ub(fix2)); % transform to unconstrained
    lgrid                   = theta_MLE2(j)-prct(j)*theta_MLE2(j);          % Smallest value grid
    ugrid                   = theta_MLE2(j)+prct(j)*theta_MLE2(j);          % Largest value grid
    grid2(j,:)              = sort([theta_MLE2(j),p(j),...                  % Set grid
                                        linspace(lgrid,ugrid,Npoints)]);    % ...
    for i = 1:Ngrid
        if j == 1
            fun             = @(x) LogLikPen([grid2(j,i);x(1:4);...         % Ro
                                                p(6:8);x0],sys);            % ...
        end
        if j == 2
            fun             = @(x) LogLikPen([x(1);grid2(j,i);x(2:4);...    % Ri
                                                p(6:8);x0],sys);            % ...
        end 
        if j == 3
            fun             = @(x) LogLikPen([x(1:2);grid2(j,i);x(3:4);...  % Cw
                                                p(6:8);x0],sys);            % ...
        end
        if j == 4
            fun             = @(x) LogLikPen([x(1:3);grid2(j,i);x(4);...    % Ci
                                                p(6:8);x0],sys);            % ...
        end
        if j == 5
            fun             = @(x) LogLikPen([x(1:4);grid2(j,i);...         % alpha
                                                p(6:8);x0],sys);            % ...
        end
        [~,fval,~,~,~,~]    = fminunc(fun,eta,options);
        PL2(j,i)            = -fval;
    end
end
%% ------------------------------------------------------------------------
% Start profile loglikelihood model M6, persistent excitation
% -------------------------------------------------------------------------
sys.u                       = u3;                                           % Find MLE
sys.y                       = y1;                                           % ...
x0                          = p(Np+Nx+Ny+1:end);                            % ...
eta                         = transform(p(fix),'LowUp',lb(fix),ub(fix));    % ...
fun                         = @(x) LogLikPen([x(1:5);p(6:8);x0],sys);       % ...
[eta_MLE,~,~,~,~,~]         = fminunc(fun,eta,options);                     % ...
theta_MLE1                  = inv_transform(eta_MLE,'LowUp',...             % ...
                                                lb(fix),ub(fix));           % ...
Npoints                     = 200;                                          % Number of grid points
prct                        = [0.1 0.1 0.1 0.1 0.5];                        % percentage of variation of the parameters
grid1                       = zeros(Np,Npoints+2);                          % Allocation grid
Ngrid                       = max(size(grid1));                             % size grid
PL1                         = zeros(Np,Ngrid);                              % Allocation profile

for j = 1:Np
    disp(['PL1: ',num2str(j),' / ',num2str(Np)])
    fix2                    = fix;                                          % Fix profile parameter
    fix2(j)                 = false;                                        % ...
    sys.fix                 = fix2;                                         % ...
    eta                     = transform(p(fix2),'LowUp',lb(fix2),ub(fix2)); % transform to unconstrained
    lgrid                   = theta_MLE1(j)-prct(j)*theta_MLE1(j);          % Smallest value grid
    ugrid                   = theta_MLE1(j)+prct(j)*theta_MLE1(j);          % Largest value grid
    grid1(j,:)              = sort([theta_MLE1(j),p(j),...                  % Set grid 
                                        linspace(lgrid,ugrid,Npoints)]);    % ...
    for i = 1:Ngrid
        if j == 1
            fun             = @(x) LogLikPen([grid1(j,i);x(1:4);...         % Ro
                                                p(6:8);x0],sys);            % ...
        end
        if j == 2
            fun             = @(x) LogLikPen([x(1);grid1(j,i);x(2:4);...    % Ri
                                                p(6:8);x0],sys);            % ...
        end
        if j == 3
            fun             = @(x) LogLikPen([x(1:2);grid1(j,i);x(3:4);...  % Cw
                                                p(6:8);x0],sys);            % ...
        end
        if j == 4
            fun             = @(x) LogLikPen([x(1:3);grid1(j,i);x(4);...    % Ci
                                                p(6:8);x0],sys);            % ...
        end
        if j == 5
            fun             = @(x) LogLikPen([x(1:4);grid1(j,i);...         % Alpha
                                                p(6:8);x0],sys);            % ...
        end
        [~,fval,~,~,~,~]    = fminunc(fun,eta,options);
        PL1(j,i)            = -fval;
    end
end

return % ------------------------------------------------------------------ PLot section below

load('M7_M6_profile.mat')
for j = 1:Np
    Xq1(j,:)                = linspace(grid1(j,1),grid1(j,end),1e4);
    Yq1(j,:)                = interp1(grid1(j,:),PL1(j,:),Xq1(j,:),'pchip');
end
for j = 1:Np
    Xq2(j,:)                = linspace(grid2(j,1),grid2(j,end),1e4);
    Yq2(j,:)                = interp1(grid2(j,:),PL2(j,:),Xq2(j,:),'pchip');
end
%% ------------------------------------------------------------------------
% Profile likelihood model M6, Figure 3.14
% -------------------------------------------------------------------------
plotheight              = 15;
plotwidth               = 15;
lw                      = 2; % LineWidth
fs                      = 12; % fontSize
fw                      = 'demi'; % fontWeight
fn                      = 'Cambria'; % fontName
dc                      = linspecer(6,'qualitative'); 
subplotsx               = 3;
subplotsy               = 2;
leftedge                = 0.6;
rightedge               = 0.6;
topedge                 = 0.5;
bottomedge              = 0.6;
spacex                  = 1.2;
spacey                  = 1.2;
sub_pos                 = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.65 0.8]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------- % 
axes('position',sub_pos{1,2});
plot(Xq1(1,:),Yq1(1,:)-max(Yq1(1,:)),'Color',dc(2,:),'Linewidth',lw);
hold on
plot(Xq2(1,:),Yq2(1,:)-max(Yq2(1,:)),'Color',dc(6,:),'Linewidth',lw);
hold on 
plot(p(1),0,'+r','Linewidth',lw+2)
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                      = get(gca,'Xtick');
nx                      = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
gc                      = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.05,0.9,'\textbf{$R_{o}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{2,2});
plot(Xq1(2,:),Yq1(2,:)-max(Yq1(2,:)),'Color',dc(2,:),'Linewidth',lw);
hold on
plot(Xq2(2,:),Yq2(2,:)-max(Yq2(2,:)),'Color',dc(6,:),'Linewidth',lw);
hold on 
plot(p(2),0,'+r','Linewidth',lw+2)
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                      = get(gca,'Xtick');
nx                      = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
gc                      = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.05,0.9,'\textbf{$R_{i}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{3,2});
plot(Xq1(3,:),Yq1(3,:)-max(Yq1(3,:)),'Color',dc(2,:),'Linewidth',lw);
hold on
plot(Xq2(3,:),Yq2(3,:)-max(Yq2(3,:)),'Color',dc(6,:),'Linewidth',lw);
hold on 
plot(p(3),0,'+r','Linewidth',lw+2)
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                      = get(gca,'Xtick');
nx                      = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
gc                      = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.05,0.9,'\textbf{$C_{w}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{1,1});
plot(Xq1(4,:),Yq1(4,:)-max(Yq1(4,:)),'Color',dc(2,:),'Linewidth',lw);
hold on
plot(Xq2(4,:),Yq2(4,:)-max(Yq2(4,:)),'Color',dc(6,:),'Linewidth',lw);
hold on 
plot(p(4),0,'+r','Linewidth',lw+2)
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                      = get(gca,'Xtick');
nx                      = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
gc                      = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.05,0.9,'\textbf{$C_{i}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{2,1});
plot(Xq1(5,:),Yq1(5,:)-max(Yq1(5,:)),'Color',dc(2,:),'Linewidth',lw);
hold on
plot(Xq2(5,:),Yq2(5,:)-max(Yq2(5,:)),'Color',dc(6,:),'Linewidth',lw);
hold on 
plot(p(5),0,'+r','Linewidth',lw+2)
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                      = get(gca,'Xtick');
nx                      = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
gc                      = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.05,0.9,'\textbf{$\alpha^{`}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')

%% ------------------------------------------------------------------------
% Plot input, Figure 3.12
% -------------------------------------------------------------------------
td                      = (0:dt:(T-1)*dt)/(3600*24);
subplotsx               = 1;
subplotsy               = 5;
leftedge                = 1;
rightedge               = 0.5;
topedge                 = 0.5;
bottomedge              = 0.9;
spacex                  = 1.2;
spacey                  = 1.1;
sub_pos                 = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.4375 0.935185185185185]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------- % 
axes('position',sub_pos{1,5});
plot(td,To,'Color','k','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
xlim([0 6])
ylabel(['\textbf{$T_{o}$}',' [$^{\circ}C$]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{1,4});
plot(td,Qo(:,1),'Color',dc(1,:),'Linewidth',lw); hold on;
plot(td,Qo(:,2),'Color',dc(3,:),'Linewidth',lw); hold on;
plot(td,Qo(:,3),'Color',dc(4,:),'Linewidth',lw); hold on;
plot(td,Qo(:,4),'Color',dc(5,:),'Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
xlim([0 6])
yy                      = get(gca,'Ytick');
ny                      = max(floor(log(abs(yy))./log(10)));
if ny >= 1
    set(gca,'YTickLabel',yy.*(10^-ny))
    text(0,1.15,['10^{',num2str(ny),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
ylabel(['\textbf{$\dot{Q}_{o_{n,e,s,w}}$}',' [W]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{1,3});
plot(td,Qo2,'Color','k','Linewidth',lw); 
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
xlim([0 6])
yy                      = get(gca,'Ytick');
ny                      = max(floor(log(abs(yy))./log(10)));
if ny >= 1
    set(gca,'YTickLabel',yy.*(10^-ny))
    text(0,1.15,['10^{',num2str(ny),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
ylabel(['\textbf{$\dot{Q}_{o}$}',' [W]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{1,2});
plot(td,Qh,'Color',dc(2,:),'Linewidth',lw); 
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
xlim([0 6])
yy                      = get(gca,'Ytick');
ny                      = max(floor(log(abs(yy))./log(10)));
if ny >= 1
    set(gca,'YTickLabel',yy.*(10^-ny))
    text(0,1.15,['10^{',num2str(ny),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
ylabel(['\textbf{$\dot{Q}_{h}$}',' [W]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{1,1});
plot(td,Qh2,'Color',dc(6,:),'Linewidth',lw); 
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
xlim([0 6])
yy                      = get(gca,'Ytick');
ny                      = max(floor(log(abs(yy))./log(10)));
if ny >= 1
    set(gca,'YTickLabel',yy.*(10^-ny))
    text(0,1.15,['10^{',num2str(ny),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
ylabel(['\textbf{$\dot{Q}_{h}$}',' [W]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
xlabel('days','fontsize',fs+2,'FontWeight',fw,'FontName',fn)
%% ------------------------------------------------------------------------
% Plot output and frequency spectrum, Figure 3.12
% -------------------------------------------------------------------------
td                      = (0:dt:(T-1)*dt)/(3600*24);
subplotsx               = 1;
subplotsy               = 2;
leftedge                = 1;
rightedge               = 0.7;
topedge                 = 0.5;
bottomedge              = 1.4;
spacex                  = 1.2;
spacey                  = 1.6;
sub_pos                 = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.4375 0.55]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------- % 
axes('position',sub_pos{1,2});
plot(td,y1,'Color',dc(2,:),'Linewidth',lw); hold on;
plot(td,y2,'Color',dc(6,:),'Linewidth',lw); 
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
xlim([0 6])
ylabel(['\textbf{$T_{i}$}',' [$^{\circ}C$]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
xlabel('days','fontsize',fs,'FontWeight',fw,'FontName',fn)
% --------------------------- % 
axes('position',sub_pos{1,1});
plot_FFT(y1,dt,dc(2,:)); hold on;
plot_FFT(y2,dt,dc(6,:)); hold on;
h1                      = 1/(12*3600);
h2                      = 1/(0.5*3600);
plot([h1 h1],[-38 65],'r','Linewidth',lw); hold on;
plot([h2 h2],[-38 30],'r','Linewidth',lw); hold on;
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
xx                      = get(gca,'Xtick');
nx                      = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.01,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
ylabel('Frequency spectrum [dB]','fontsize',fs,'FontWeight',fw,'FontName',fn)
xlabel('Frequency [Hz]','fontsize',fs,'FontWeight',fw,'FontName',fn)
text(h1,70,'12h','Color','r','FontSize',fs,'FontWeight',fw,'FontName',fn,'HorizontalAlignment','center')
text(h2,35,'0.5h','Color','r','FontSize',fs,'FontWeight',fw,'FontName',fn,'HorizontalAlignment','center')
