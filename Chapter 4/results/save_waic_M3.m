% Compute model M3

clear all
close all
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');

%% ------------------------------------------------------------------------
% Import data for model M3
% -------------------------------------------------------------------------
load('M3_5alpha_chains.mat')
impM3                       = load('M3_5alpha_chains.mat');
stop                        = sys.MC;
start                       = 501;                                          % burnin
chain1                      = reshape(p_ALL(1,:,:),[Nukwn,stop]);
chain2                      = reshape(p_ALL(2,:,:),[Nukwn,stop]);
chain3                      = reshape(p_ALL(3,:,:),[Nukwn,stop]);
chain5                      = reshape(p_ALL(5,:,:),[Nukwn,stop]);
chain6                      = reshape(p_ALL(6,:,:),[Nukwn,stop]);

%% ------------------------------------------------------------------------
% Inverse transform Markov chains to constrained space
% -------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------
% Model configuration
% -------------------------------------------------------------------------
sys.fun                     = @M3_5alpha;                                   % Model M3
thetaAll_M3                 = [theta1_M3(:,start:stop) ...                  % gather all Markov chains
                               theta2_M3(:,start:stop) ...                  % ...
                               theta3_M3(:,start:stop) ...                  % ...
                               theta4_M3(:,start:stop) ...                  % ...
                               theta5_M3(:,start:stop)];                    % ...
[Utheta,~,IC]               = unique(thetaAll_M3','rows');                  % Take only unique parameter values
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
    par(~fix)               = y(1);                                         % Pass fixed parameters
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
