% SImulate from posterior distribution, model M3
% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function

clear all
close all
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');
addpath(userpath, 'export_fig');

%% ------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------
% Inverse transform Markov chains to constrained space
% -------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------
% Plot configuration
% -------------------------------------------------------------------------
lw                          = 2;                                            % LineWidth
fs                          = 14;                                           % fontSize
fw                          = 'demi';                                       % fontWeight
fn                          = 'Times New Roman';                            % fontName
dcolors                     = linspecer(9,'qualitative');                   % Distinguishable colors

%% ------------------------------------------------------------------------
% Model configuration
% -------------------------------------------------------------------------
start                       = 501;                                          % burnin
sys.fun                     = @M4_5alpha;                                   % Model M4
thetaAll_M4                 = [theta1_M4(:,start:stop) ...                  % gather all Markov chains
                               theta2_M4(:,start:stop) ...                  % ...
                               theta3_M4(:,start:stop) ...                  % ...
                               theta4_M4(:,start:stop) ...                  % ...
                               theta5_M4(:,start:stop)];                    % ...
[Utheta,~,IC]               = unique(thetaAll_M4','rows');                  % Take only unique parameter values
Utheta                      = Utheta';                                      % ...
count                       = accumarray(IC, 1);                            % cumulated number of unique values
cutp                        = 5001:6441;                                    % Index measure
cutu                        = 5000:6440;                                    % Index input
Tp                          = length(cutu);                                 % length input
ypt                         = Tzone1_all(cutp);                             % measured output
up                          = [Text(cutu) Tzone2_all(cutu) Sr_west(cutu)... % Input vector
                               Sr_east(cutu) Sr_south(cutu) Qhvac1(cutu)]'; % ...
afoh                        = zeros(Nu,Tp-1);                               % By default zero order hold 'zoh' (4.8)
if strcmp(sys.ahold,'foh')                                                  % if first order hold
    afoh                    = (up(:,2:end) - up(:,1:end-1))./dt;            % (4.8)
end
N                           = max(size(Utheta));                            % number of unique values
Nall                        = max(size(thetaAll_M4));                       % total size of sample
yp                          = zeros(Nall,Tp-1);                             % Allocation for simulated posterior 
nn                          = 1;                                            % init count

for n = 1:N
    par                     = zeros(Np+2*Nx+Ny,1);                          % Init parameter vector
    par(~fix)               = [2.52e-3 ; 2.52e-4 ; 1e-6 ; y(1)];            % Pass fixed parameters
    par(fix)                = Utheta(:,n);                                  % Pass parameters
    [~,xfit]                = Residuals(par,sys);                           % get last state estimates from the fitting data set
    xp                      = xfit(:,end);                                  % ...
    model                   = feval(sys.fun,par,sys);                       % Model M4 evaluation
    A                       = model.A;                                      % state matrix
    B                       = model.B;                                      % input matrix
    F                       = expm(A*dt);                                   % Discretization (4.7)
    G                       = zeros(Nx,2*Nu);                               % ...
    bis                     = F-eye(Nx);                                    % ...
    G(:,1:Nu)               = A\bis*B;                                      % ...
    if strcmp(sys.ahold,'foh')                                              % ...
        G(:,Nu+1:2*Nu)      = A\(-A\bis + F*dt)*B;                          % ...
    end                                                                     % ...
    Qstd                    = diag(par(Np+1:Np+Nx));                        % get std process noise
    Rstd                    = par(Np+Nx+1:Np+Nx+Ny);                        % get std measurement noise
    for kk = 1:count(n)                                                     % For count(n) similar parameters, we simulate count(n)
        for k = 1:(Tp-1) % afoh [Nu x Tp-1]
            xp              = F*xp + G(:,1:Nu)*(afoh(:,k)*dt + up(:,k)) ... % Simulation
                              - G(:,Nu+1:2*Nu)*afoh(:,k) + Qstd*randn(Nx,1);% ...
            yp(nn,k)        = sys.H*xp + Rstd.*randn(Ny,1);                 % ...
        end
        nn                  = nn + 1;                                       % keep count
    end
    disp([num2str(n), '/',num2str(N)])                                      % display iterations
end

MinMax                      = zeros(2,Tp-1);                                % Allocation for min/max simulated output
for k = 1:(Tp-1)                                                            % At each time step, find min and max
    MinMax(1,k)             = max(yp(:,k));                                 % ...
    MinMax(2,k)             = min(yp(:,k));                                 % ...
end

%% ------------------------------------------------------------------------
% Time in MATLAB date format
% -------------------------------------------------------------------------
startDate                   = datenum('13-05-2014 17:30:00', 'dd-mm-yy HH:MM:SS');
endDate                     = datenum('23-05-2014 17:30:00', 'dd-mm-yy HH:MM:SS'); 
xData                       = linspace(startDate,endDate,Tp-1);

%% ------------------------------------------------------------------------
% PLot simulation from posterior
% -------------------------------------------------------------------------
figure('visible','on','units','normalized','outerposition',[0 0.0648148148148148 0.4375 0.5]);
X                           = [xData,fliplr(xData)];
Y                           = [MinMax(2,:),fliplr(MinMax(1,:))];
h                           = fill(X,Y,dcolors(2,:),'edgecolor','none');
set(h,'facealpha',0.8)
hold on
plot(xData,ypt(1:Tp-1),'r','linewidth',1.5)
datetick('x','dd','keepticks')
axis tight
set(gca,'Box','off','FontWeight',fw,'FontSize',fs,'FontName',fn,'Linewidth',lw)
set(gcf,'color','w');
ylabel('\textbf{$^{\circ}C$}','FontSize',16,'FontWeight',fw,'FontName',fn,'Units','Normalized','interpreter','latex');
text(1.02,0.02,'May','FontSize',fs,'FontWeight',fw,'Units','Normalized','FontName',fn);

% save('simM4.mat','X','Y','yp','ypt','xData')