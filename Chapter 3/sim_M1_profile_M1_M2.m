% Plot
%   - Figure 3.7: Inputs and Outputs data
%   - Figure 3.8: profile M1 with data simulated with M1
%   - Figure 3.9: surface profile Ro Rwo with data simulated with M1 
%   - Figure 3.10: profile M2 with data simulated with M1

% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig

clear all;
close all;
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');
%% ------------------------------------------------------------------------
% Import data | Inputs and Outputs definition
% Import Twin Houses experiment data available at 
% https://pure.strath.ac.uk/portal/en/datasets/twin-houses-empirical-validation-dataset-experiment-2(94559779-e781-4318-8842-80a2b1201668).html
% The data of the Twin Houses experiment are used for the simulated examples
% -------------------------------------------------------------------------
dt                          = 10*60;                                        % sampling time [s]
imp                         = importdata('Twin_house_exp2_10min.xlsx');     % Import data
cut                         = 1:6400;                                       % part of data
T                           = length(cut);                                  % Length of data
Tliving010                  = imp.data(:,6);                                % Temperature living room at 10 cm from the ground [°C]
Tliving110                  = imp.data(:,7);                                % Temperature living room at 110 cm from the ground [°C]
Tliving170                  = imp.data(:,10);                               % Temperature living room at 170 cm from the ground [°C]
Tliving                     = mean([Tliving010,Tliving110,Tliving170],2);   % Mean temperature of the living room  [°C]
Tcorridor                   = imp.data(:,11);                               % Temperature corridor [°C]
Tbath010                    = imp.data(:,12);                               % Temperature in the bathroom at 10 cm from the ground [°C]
Tbath110                    = imp.data(:,13);                               % Temperature in the bathroom at 110 cm from the ground [°C]
Tbath170                    = imp.data(:,14);                               % Temperature in the bathroom at 170 cm from the ground [°C]
Tbath                       = mean([Tbath010,Tbath110,Tbath170],2);         % Mean temperature of the bathroom [°C]
Tchild010                   = imp.data(:,15);                               % Temperature in the child room at 10 cm from the ground [°C]
Tchild110                   = imp.data(:,16);                               % Temperature in the child room at 110 cm from the ground [°C]
Tchild170                   = imp.data(:,17);                               % Temperature in the child room at 170 cm from the ground [°C]
Tchild                      = mean([Tchild010,Tchild110,Tchild170],2);      % Mean temperaure of the child room [°C]
Text                        = imp.data(:,54);                               % Outdoor air temperature [°C]
Sgliving                    = 33.65;                                        % ground surface living room 
Sgchild                     = 11.19;                                        % ground surface child room
Sgbath                      = 6.92;                                         % ground surface bathroom
Sgcorridor                  = 5.46;                                         % ground surface corridor
Sg1                         = Sgliving+ Sgchild + Sgbath + Sgcorridor;      % ground surface south zone
Vspaces                     = [Sgliving;Sgchild;Sgbath;Sgcorridor]./Sg1;    % Normalized surfaces of the south zone
Tzone1_all                  = [Tliving Tbath Tchild Tcorridor]*Vspaces;     % south zone temperature, weighted mean according to ground surfaces [°C]
To                          = Text(cut);                                    % Take only a part of the data
Ti                          = Tzone1_all(cut);                              % Take only a part of the data
u                           = [To Ti]';                                     % Input vector
%% ------------------------------------------------------------------------
% Simulate known model 
% -------------------------------------------------------------------------
Np                          = 5;                                            % Number of R,C parameters
Nx                          = 1;                                            % Number of states 
Ny                          = 1;                                            % Number of outputs
Nu                          = 2;                                            % Number of inputs
fix                         = [true(Np,1) ; false(2*Nx+Ny,1)];              % logical indexing of fixed parameters
sys.dim                     = [Np Nx Ny Nu];                                % vector of dimensions
sys.fix                     = fix;                                          % Pass logical indexing into the structure
sys.fun                     = @M1;                                          % Simulation with model M1
sys.Psqrt                   = 0.5;                                          % Initial std of state                                    
sys.sc                      = 1e8;                                          % scaling of thermal capacity
sys.dt                      = dt;                                           % sampling time

S                           = 10;                                           % surface [m²], hard coding in @M1
e                           = 0.2;                                          % width [m]
Rso                         = 1/(25*S);                                     % Outdoor convective resistance [K/W] (3.39)
Rsi                         = 1/(8*S);                                      % Indoor convective resistance [K/W] (3.39)
lambda                      = 2;                                            % Thermal conductivity [W/(m.K)]
rho                         = 2400;                                         % density [kg.m^-3]
c                           = 1000;                                         % specific heat capacity [J/(Kg.K)]
Rw                          = e/(S*lambda);                                 % thermal resistance of the concrete wall [K/W] (3.36)
Rwo                         = Rw/2;                                         % symmetric wall (3.37)                                         
Rwi                         = Rw/2;                                         % ...
Cw                          = rho*S*e*c;                                    % thermal capacity of the concrete wall [K/W] (3.38)
Ro                          = Rso+Rwo;                                      % thermal resistance from the outside to the middle of the wall [K/W]
Ri                          = Rsi+Rwi;                                      % thermal resistance from the middle of the wall to the inside [K/W]
sigw                        = 1e-2;                                         % std of process noise
sigv                        = 1e-2;                                         % std of measurement noise
xw0                         = 15;                                           % initial state

% M1                            lb     p     ub                            
gather(:,1)                 = [1e-4 ; Rso  ; 1];
gather(:,2)                 = [1e-4 ; Rwo  ; 1];
gather(:,3)                 = [1e-4 ; Rwi  ; 1];
gather(:,4)                 = [1e-4 ; Rsi  ; 1];
gather(:,5)                 = [1e-3 ; Cw/sys.sc ; 1];
gather(:,6)                 = [1e-8 ; sigw ; 1];
gather(:,7)                 = [1e-8 ; sigv ; 1];
gather(:,8)                 = [10   ; xw0  ; 40];
p                           = gather(2,:)';                                 % parameter vector of M1
lb                          = gather(1,:)';                                 % lower bounds
ub                          = gather(3,:)';                                 % upper bounds
sys.lb                      = lb;                                           % pass lb to the structure
sys.ub                      = ub;                                           % pass ub to the structure

x                           = p(Np+Nx+Ny+1:end);                            % get parameters for the simulation 
Qstd                        = p(Np+1:Np+Nx);                                % ...
Rstd                        = p(Np+Nx+1:Np+Nx+Ny);                          % ...
y                           = zeros(Ny,T);                                  % Allocation for the simulated output
xS                          = zeros(Nx,T);                                  % Allocation for the simulated state
model                       = feval(sys.fun,p(1:Np),sys);                   % @M1
A                           = model.A;                                      % state matrix
B                           = model.B;                                      % input matrix
C                           = model.C;                                      % output matrix
D                           = model.D;                                      % feedthrough matrix
F                           = expm(A*dt);                                   % discretization
G                           = A\(F-eye(Nx))*B;                              % ...
for k = 1:T                                                                 % Simulate model M1
    xS(:,k)                 = x;                                            % ...
    y(k)                    = C*x + D*u(:,k) + Rstd*randn(Ny,1);            % ...
    x                       = F*x + G*u(:,k) + Qstd*randn(Nx,1);            % ...
end
cut2                        = 2890:4829;                                    % subset of the simulated data
sys.u                       = u(:,cut2);                                    % Pass inputs to the structure
sys.y                       = y(cut2);                                      % Pass output to the structure
x0                          = xS(:,cut2(1));                                % Initial state of the subset

%% ------------------------------------------------------------------------
% Find MLE of model M1
% -------------------------------------------------------------------------
options                     = optimoptions(@fminunc,...                     % Quasi Newton with central difference and BFGS
                                    'Display','none',...
                                    'Algorithm','quasi-newton',...
                                    'FiniteDifferenceType','central',...
                                    'HessUpdate','bfgs',...
                                    'MaxIterations',1e3,...
                                    'MaxFunctionEvaluations',1e3,...
                                    'OptimalityTolerance',1e-9,...
                                    'UseParallel',true,...
                                    'StepTolerance',1e-9);
eta                         = transform(p(fix),'LowUp',lb(fix),ub(fix));    % transform to unconstrained space
fun                         = @(x) LogLikPen([x(1:5);p(6:7);x0],sys);       % Penalyzed LogLikelihood without gradient 
[eta_MLE1,~,~,~,~,~]        = fminunc(fun,eta,options);                     % optimize
theta_MLE1                  = inv_transform(eta_MLE1,'LowUp', ...           % inverse transform to theta
                                                lb(fix),ub(fix));           % ...

%% ------------------------------------------------------------------------
% Start: Profile likelihood model M1
% -------------------------------------------------------------------------
Npoints                     = 100;                                          % Number of points for the profiles
prct                        = [0.1 0.1 0.1 0.1 0.03];                       % +/- of theta_MLE
grid1                       = zeros(Np,Npoints+2);                          % Allocation grid values
Ngrid                       = max(size(grid1));                             % length grid
PL1                         = zeros(Np,Ngrid);                              % Allocation profile loglikelihood

for j = 1:Np
    disp(['PL1: ',num2str(j),' / ',num2str(Np)])                            % display iteration
    fix2                    = fix;                                          % for each j, fix the corresponding parameter to a grid value
    fix2(j)                 = false;                                        % ...
    sys.fix                 = fix2;                                         % ...
    eta                     = transform(p(fix2),'LowUp',lb(fix2),ub(fix2)); % transform to unconstrained space
    lgrid                   = theta_MLE1(j)-prct(j)*theta_MLE1(j);          % Smallest value grid
    ugrid                   = theta_MLE1(j)+prct(j)*theta_MLE1(j);          % largest value grid
    grid1(j,:)              = sort([theta_MLE1(j),p(j),...                  % create grid
                                        linspace(lgrid,ugrid,Npoints)]);    % ...
    for i = 1:Ngrid
        if j == 1
            fun             = @(x) LogLikPen([grid1(j,i);x(1:4);...         % Rso
                                                p(6:7);x0],sys);            % ...
        end
        if j == 2
            fun             = @(x) LogLikPen([x(1);grid1(j,i);x(2:4);...    % Rwo
                                                p(6:7);x0],sys);            % ...
        end
        if j == 3
            fun             = @(x) LogLikPen([x(1:2);grid1(j,i);x(3:4);...  % Rwi
                                                p(6:7);x0],sys);            % ...
        end
        if j == 4
            fun             = @(x) LogLikPen([x(1:3);grid1(j,i);x(4);...    % Rsi
                                                p(6:7);x0],sys);            % ...
        end
        if j == 5
            fun             = @(x) LogLikPen([x(1:4);grid1(j,i);...         % Cw
                                                p(6:7);x0],sys);            % ...
        end
        [~,fval,~,~,~,~]    = fminunc(fun,eta,options);                     % Unconstrained optimization 
        PL1(j,i)            = -fval;                                        % Save profile values
    end
end
%% ------------------------------------------------------------------------
% Find MLE of model M2
% -------------------------------------------------------------------------
Np                          = 3;                                            % Number of R,C parameters
Nx                          = 1;                                            % Number of states
Ny                          = 1;                                            % Number of outputs
Nu                          = 2;                                            % Number of inputs
fix                         = [true(Np,1) ; false(2*Nx+Ny,1)];              % logical indexing of fixed parameters
sys.dim                     = [Np Nx Ny Nu];                                % vector of dimensions
sys.fix                     = fix;                                          % Pass to the structure
sys.fun                     = @M2;                                          % Model M2
sys.Psqrt                   = 0.5;                                          % Initial state std
%                               lb     p     ub
gather2(:,1)                = [1e-4 ; Rso+Rwo  ; 1];
gather2(:,2)                = [1e-4 ; Rsi+Rwi  ; 1];
gather2(:,3)                = [1e-3 ; Cw/sys.sc ; 1];
gather2(:,4)                = [1e-8 ; sigw ; 1];
gather2(:,5)                = [1e-8 ; sigv ; 1];
gather2(:,6)                = [1   ; xw0  ; 40];
p                           = gather2(2,:)';                                % parameters of M2
lb                          = gather2(1,:)';                                % lower bounds
ub                          = gather2(3,:)';                                % upper bounds
sys.lb                      = lb;                                           % Pass to the structure
sys.ub                      = ub;                                           % Pass to the structure

eta                         = transform(p(fix),'LowUp',lb(fix),ub(fix));    % Convert to unconstrained space
fun                         = @(x) LogLikPen([x(1:3);p(4:5);x0],sys);       % Penalyzed LogLikelihood without gradient
[eta_MLE,~,~,~,~,~]         = fminunc(fun,eta,options);                     % Unconstrained optimization 
theta_MLE2                  = inv_transform(eta_MLE,'LowUp',...             % inverse transform to theta_MLE(M2)
                                                lb(fix),ub(fix));           % ...
%% ------------------------------------------------------------------------
% Start profile loglikelihood model M2
% -------------------------------------------------------------------------
Npoints                     = 100;                                          % Number of points for profile
prct                        = [0.01 0.01 0.03];                             % percentage of theta_MLE(M2) 
grid2                       = zeros(Np,Npoints+2);                          % Allocation grid 
Ngrid                       = max(size(grid2));                             % length of grid
PL2                         = zeros(Np,Ngrid);                              % Allocation profile M2

for j = 1:Np
    disp(['PL2: ',num2str(j),' / ',num2str(Np)])                            % display iteration
    fix2                    = fix;                                          % for each j, fix the corresponding parameter to a grid value
    fix2(j)                 = false;                                        % ...
    sys.fix                 = fix2;                                         % ...
    eta                     = transform(p(fix2),'LowUp',lb(fix2),ub(fix2)); % transform to unconstrained parameters
    lgrid                   = theta_MLE2(j)-prct(j)*theta_MLE2(j);          % Smallest grid value  
    ugrid                   = theta_MLE2(j)+prct(j)*theta_MLE2(j);          % largest grid value 
    grid2(j,:)              = sort([theta_MLE2(j),p(j),...                  % create grid
                                    linspace(lgrid,ugrid,Npoints)]);        % ...
    for i = 1:Ngrid
        if j == 1
            fun             = @(x) LogLikPen([grid2(j,i);x(1:2);...         % Ro
                                                p(4:5);x0],sys);            % ...
        end
        if j == 2
            fun             = @(x) LogLikPen([x(1);grid2(j,i);x(2);...      % Ri
                                                p(4:5);x0],sys);            % ...
        end
        if j == 3
            fun             = @(x) LogLikPen([x(1:2);grid2(j,i);...         % Cw
                                                p(4:5);x0],sys);            % ...
        end
        [~,fval,~,~,~,~]    = fminunc(fun,eta,options);
        PL2(j,i)            = -fval;                                        % Save profile values
    end
end
%% ------------------------------------------------------------------------
% Start surface profile, model M1, Rso and Rwo 
% -------------------------------------------------------------------------
Np                          = 5;                                            % Number of R,C parameters
Nx                          = 1;                                            % Number of states
Ny                          = 1;                                            % Number of outputs
Nu                          = 2;                                            % Number of inputs
fix                         = [true(Np,1) ; false(2*Nx+Ny,1)];              % Logical indexing of fixed parameters
sys.dim                     = [Np Nx Ny Nu];                                % vector of dimensions
sys.fix                     = fix;                                          % Passing to the structure
sys.fun                     = @M1;                                          % Model M1
sys.Psqrt                   = 0.5;                                          % std of initial state
sys.sc                      = 1e8;                                          % scaling of thermal capacities
sys.dt                      = dt;                                           % sampling time
%                               lb     p     ub
gather(:,1)                 = [1e-4 ; Rso  ; 1];
gather(:,2)                 = [1e-4 ; Rwo  ; 1];
gather(:,3)                 = [1e-4 ; Rwi  ; 1];
gather(:,4)                 = [1e-4 ; Rsi  ; 1];
gather(:,5)                 = [1e-3 ; Cw/sys.sc ; 1];
gather(:,6)                 = [1e-8 ; sigw ; 1];
gather(:,7)                 = [1e-8 ; sigv ; 1];
gather(:,8)                 = [10   ; xw0  ; 40];   
p                           = gather(2,:)';                                 % parameters
lb                          = gather(1,:)';                                 % lower bounds
ub                          = gather(3,:)';                                 % upper bounds
sys.lb                      = lb;                                           % Passing to the structure
sys.ub                      = ub;                                           % ...

Npoints                     = 100;                                          % Number of point for surface loglikelihood
prct                        = 0.02;                                         % percentage of variation 
lgrid                       = theta_MLE1(1:2)-prct.*theta_MLE1(1:2);        % Smallest value grid
ugrid                       = theta_MLE1(1:2)+prct.*theta_MLE1(1:2);        % Largest value grid
grid11                      = sort([theta_MLE1(1),p(1),linspace(...         % grid Rso
                                        lgrid(1),ugrid(1),Npoints)]);       % ...
grid22                      = sort([theta_MLE1(2),p(2),linspace(...         % grid Rwo
                                        lgrid(2),ugrid(2),Npoints)]);       % ...
Ngrid                       = max(size(grid11));                            % length grid
PL3                         = zeros(Ngrid,Ngrid);                           % Allocation surface values

fix2                        = fix;                                          % fix Rso and Rwo
fix2(1:2)                   = false;                                        % ...
sys.fix                     = fix2;                                         % Passing to the structure
eta                         = transform(p(fix2),'LowUp',lb(fix2),ub(fix2)); % transform to unconstrained parameters
k                           = 1;                                            % init count
for i = 1:Ngrid                                                             % Surface loglikelihood
    for j = 1:Ngrid                                                         % ...
        disp([num2str(k),' / ',num2str(Ngrid^2)])                           % ...
        fun                 = @(x) LogLikPen([grid11(i);grid22(j);...       % ...
                                                x(1:3);p(6:7);x0],sys);     % ...
        [~,fval,~,~,~,~]    = fminunc(fun,eta,options);                     % ...
        PL3(i,j)            = -fval;                                        % ...
        k                   = k + 1;                                        % ...
    end                                                                     % ...
end                                                                         % ...

return % ------------------------------------------------------------------ PLot section below

load('M1_M2_profile_V2.mat')
%% ------------------------------------------------------------------------
% Start: Surface plot, Figure 3.9
% -------------------------------------------------------------------------
plotheight                  = 15;
plotwidth                   = 15;
lw                          = 2;                                            % LineWidth
fs                          = 12;                                           % fontSize
fw                          = 'demi';                                       % fontWeight
fn                          = 'Cambria';                                    % fontName
figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.4375 0.6]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
[X,Y]                       = meshgrid(grid11,grid22);
gridi1                      = linspace(grid11(1),grid11(end),1e3);
gridi2                      = linspace(grid22(1),grid22(end),1e3);
[Xi,Yi]                     = meshgrid(gridi1,gridi2);
PLi                         = interp2(X,Y,PL3,Xi,Yi,'linear');
Z                           = PLi-max(PLi);
colormap(jet) 
surf(Xi,Yi,Z,'EdgeColor','interp','FaceColor','interp')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
colorbar
view(0,90)
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1,-0.026,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
yy                          = get(gca,'Ytick');
ny                          = max(floor(log(abs(yy))./log(10)));
ny                          = ny - 1; 
if ny <= -1
    set(gca,'YTickLabel',yy.*(10^-ny))
    text(-0.005,1.05,['10^{',num2str(ny),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
xlabel('\textbf{$R_{so}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
ylabel('\textbf{$R_{wo}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')

%% ------------------------------------------------------------------------
% Profile M1, Figure 3.8 
% -------------------------------------------------------------------------
subplotsx                   = 3;
subplotsy                   = 2;
leftedge                    = 0.8;
rightedge                   = 0.8;
topedge                     = 0.5;
bottomedge                  = 0.6;
spacex                      = 1.3;
spacey                      = 1.2;
sub_pos                     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.65 0.8]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------- % 
axes('position',sub_pos{1,2});
plot(grid1(1,:),PL1(1,:)-max(PL1(1,:)),'Color','b','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
plot(Rso,0,'+r','Linewidth',lw+2)
hold on 
gc                          = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.5,0.65,'\textbf{$R_{so}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{2,2});
plot(grid1(2,:),PL1(2,:)-max(PL1(2,:)),'Color','b','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
plot(Rwo,0,'+r','Linewidth',lw+2)
hold on 
gc                          = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.5,0.65,'\textbf{$R_{wo}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{3,2});
plot(grid1(3,:),PL1(3,:)-max(PL1(3,:)),'Color','b','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
plot(Rwi,0,'+r','Linewidth',lw+2)
hold on 
gc                          = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.5,0.65,'\textbf{$R_{wi}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{1,1});
plot(grid1(4,:),PL1(4,:)-max(PL1(4,:)),'Color','b','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
plot(Rsi,0,'+r','Linewidth',lw+2)
hold on 
gc                          = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.5,0.65,'\textbf{$R_{si}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{2,1});
plot(grid1(5,:),PL1(5,:)-max(PL1(5,:)),'Color','b','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
plot(Cw/sys.sc,0,'+r','Linewidth',lw+2)
hold on 
gc = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.45,0.65,'\textbf{$10^8 C_{w}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')

%% ------------------------------------------------------------------------
% Profile M2, Figure 3.10 
% -------------------------------------------------------------------------
subplotsx                   = 3;
subplotsy                   = 1;
leftedge                    = 0.6;
rightedge                   = 0.6;
topedge                     = 0.5;
bottomedge                  = 1;
spacex                      = 1;
spacey                      = 1.2;
sub_pos                     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.65 0.5]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------- % 
axes('position',sub_pos{1,1});
plot(grid2(1,:),PL2(1,:)-max(PL2(1,:)),'Color','b','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
plot(Ro,0,'+r','Linewidth',lw+2)
hold on 
gc                          = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.5,0.65,'\textbf{$R_{o}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{2,1});
plot(grid2(2,:),PL2(2,:)-max(PL2(2,:)),'Color','b','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
plot(Ri,0,'+r','Linewidth',lw+2)
hold on 
gc                          = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.5,0.65,'\textbf{$R_{i}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{3,1});
plot(grid2(3,:),PL2(3,:)-max(PL2(3,:)),'Color','b','Linewidth',lw);
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylim([-3 0])
xx                          = get(gca,'Xtick');
nx                          = max(floor(log(abs(xx))./log(10)));
if nx <= -1
    set(gca,'XTickLabel',xx.*(10^-nx))
    text(1.03,0.025,['10^{',num2str(nx),'}'],'FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);
end
hold on 
plot(Cw/sys.sc,0,'+r','Linewidth',lw+2)
hold on 
gc                          = get(gca);
plot([gc.XLim(1) gc.XLim(2)],[-1.92 -1.92],'--k','Linewidth',lw)
text(0.5,0.65,'\textbf{$C_{w}$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex')

%% ------------------------------------------------------------------------
% Plot Input/output Figure 3.7 
% -------------------------------------------------------------------------
td                          = (0:dt:(T2-1)*dt)/(3600*24);
subplotsx                   = 1;
subplotsy                   = 3;
leftedge                    = 1.2;
rightedge                   = 0.7;
topedge                     = 0.5;
bottomedge                  = 1.4;
spacex                      = 1.2;
spacey                      = 1.6;
sub_pos                     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.4375 0.55]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------- % 
axes('position',sub_pos{1,3});
plot(td,y2,'Color','k','Linewidth',lw); 
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylabel(['\textbf{$\dot{Q}_{i}$}',' [W]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{1,2});
plot(td,u2(1,:),'Color','k','Linewidth',lw); 
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylabel(['\textbf{$T_{o}$}',' [$^{\circ}C$]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
% --------------------------- % 
axes('position',sub_pos{1,1});
plot(td,u2(2,:),'Color','k','Linewidth',lw); 
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
axis tight
ylabel(['\textbf{$T_{i}$}',' [$^{\circ}C$]'],'fontsize',fs+2,'FontWeight',fw,'FontName',fn,'interpreter','latex')
xlabel('days','fontsize',fs,'FontWeight',fw,'FontName',fn)
