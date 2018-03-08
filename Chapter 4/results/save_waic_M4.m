% Compute model M4

clear all
close all
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');

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
cut                         = 2890:5000;                                    % identification data set
cutp                        = 5000:6441;                                    % validation data set
y                           = Tzone1_all(cut);                              % identification output
yp                          = Tzone1_all(cutp);                             % validation output
u                           = [Text(cut) Tzone2_all(cut) Sr_west(cut) ...   % Identification input vector
                               Sr_east(cut) Sr_south(cut) Qhvac1(cut)]';    % ...
up                          = [Text(cutp) Tzone2_all(cutp) Sr_west(cutp) ...% Validation input vector
                               Sr_east(cutp) Sr_south(cutp) Qhvac1(cutp)]'; % ...
T                           = length(cut);                                  % Data length identification set
Tp                          = length(cutp);                                 % Data length validation set
N                           = max(size(Utheta));                            % number of unique values
lppdp                       = zeros(N,Tp);                                  % Allocation for loglikelihood validation

%% ------------------------------------------------------------------------
% Compute WAIC
% -------------------------------------------------------------------------
for n = 1:N
    par                     = zeros(Np+2*Nx+Ny,1);                          % Init parameter vector
    par(~fix)               = [2.52e-3 ; 2.52e-4 ; 1e-6 ; y(1)];            % Pass fixed parameters    
    par(fix)                = Utheta(:,n);                                  % Pass parameters
    sys.y                   = y;                                            % get last state estimates from the fitting data set
    sys.u                   = u;                                            % ...
    [~,xfit,~]              = Residuals(par,sys);                           % ...
    par(end-Nx+1:end)       = xfit(:,end);                                  % Compute loglikelihood validation data set
    sys.y                   = yp;                                           % ...
    sys.u                   = up;                                           % ...
    [~,~,lppdp(n,:)]        = Residuals(par,sys);                           % ...
    disp([num2str(n),' / ',num2str(N)])                                     % display iterations
end

lppdAllp                    = repelem(lppdp,count,1);                       % Take into account similar values
mPost                       = sum(mean(lppdAllp,1));                        % mean 
vPost                       = sum(var(lppdAllp,1));                         % variance
WAIC                        = -2*(mPost - vPost);                           % WAIC
