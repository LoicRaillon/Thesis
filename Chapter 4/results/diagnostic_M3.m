% Compute:
%   - acceptance ratio
%   - Potential Scale Reduction Factor, R
%   - Effective sample size, neff
%   - Integrated autocorrelation time, IACT
%   - Time cosntants
%   - 95% probability mass of posterior distribution 

% Functions used: 
% - M3_5alpha: Model M3
% - inv_transform: theta = f^-1(eta) 
% - IACT_plot: Compute IACT and plot autocorrelation of the Markov Chains
% - psrf: Gelman-Rubin Potential Scale Reduction Factor, fonction accessible in the MCMC Diagnostics for Matlab at: http://becs.aalto.fi/en/research/bayes/mcmcdiag/  

clear all
close all
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');
addpath(userpath, 'export_fig');

%% ------------------------------------------------------------------------
% Import data for model M3
% -------------------------------------------------------------------------
load('M3_5alpha_chains.mat')
impM3                       = load('M3_5alpha_chains.mat');
stop                        = sys.MC;
start                       = 501; % Burn-in 
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
% Average acceptance ratio of the 5 Markov chains
% -------------------------------------------------------------------------
ratio                       = mean(sum(acc_ALL(:,start:stop),2)./length(acc_ALL(:,start:stop)));

%% ------------------------------------------------------------------------
% Convergence diagnosis, section 4.2.2.5
% -------------------------------------------------------------------------
% R : Potential scale factor reduction
% neff : Effective sample size

% First split each chains in two to have 2M chains
% Remove burn-in
% M chains
% N length chains
% Np number of parameters
% X = [NxNpxM] matrix 
% N = [MC - burnin]/2
N                           = (sys.MC - start+1)/2;
ch1                         = theta1_M3(:,start+1:start+N)';
ch2                         = theta1_M3(:,start+N:end)';
ch3                         = theta2_M3(:,start+1:start+N)';
ch4                         = theta2_M3(:,start+N:end)';
ch5                         = theta3_M3(:,start+1:start+N)';
ch6                         = theta3_M3(:,start+N:end)';
ch7                         = theta4_M3(:,start+1:start+N)';
ch8                         = theta4_M3(:,start+N:end)';
ch9                         = theta5_M3(:,start+1:start+N)';
ch10                        = theta5_M3(:,start+N:end)';
X                           = cat(3,ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9,ch10);
[R,neff,~,~,~,~,~]          = psrf(X);

%% ------------------------------------------------------------------------
% Integrated autocorrelation time (IACT), section 4.2.2.4
% -------------------------------------------------------------------------
nlags                   = 144; % Number of lags
ncol                    = 3; % Number of columns in the figure
iact1                   = IACT_plot(theta1_M3,start,nlags,ncol);
iact2                   = IACT_plot(theta2_M3,start,nlags,ncol);
iact3                   = IACT_plot(theta3_M3,start,nlags,ncol);
iact4                   = IACT_plot(theta4_M3,start,nlags,ncol);
iact5                   = IACT_plot(theta5_M3,start,nlags,ncol);
iact                    = [iact1 iact2 iact3 iact4 iact5];
mmiact                  = [min(iact,[],2) max(iact,[],2)]; % Min/Max IACT

%% ------------------------------------------------------------------------
% Time constants of model M3
% -------------------------------------------------------------------------
thetaAll_M3             = [theta1_M3(:,start:stop) theta2_M3(:,start:stop) theta3_M3(:,start:stop) theta4_M3(:,start:stop) theta5_M3(:,start:stop)];
Utheta                  = unique(thetaAll_M3','rows','stable')'; % Remove similar theta
nUni                    = max(size(Utheta));
tau                     = zeros(Nx,nUni);
for i = 1:nUni
    ssm                 = feval(@M3_5alpha,Utheta(:,i),sys);
    tau(:,i)            = (-1./eig(ssm.A))./3600;
end
% Min/Max time constants
[min(tau(1,:)) max(tau(1,:))]
[min(tau(2,:)) max(tau(2,:))]
[min(tau(3,:)) max(tau(3,:))] 

%% ------------------------------------------------------------------------
% 95% of the posterior probability distribution
% -------------------------------------------------------------------------   
ciAll                   = zeros(Nukwn,2);
for i = 1:Nukwn
    [f,x]               = ecdf(thetaAll_M3(i,:));
    [~,idxL]            = min(abs(f-0.025));
    [~,idxU]            = min(abs(f-0.975));
    ciAll(i,:)          = [x(idxL) x(idxU)];
end
[f,x]                   = ecdf(thetaAll_M3(1,:)+thetaAll_M3(2,:));
[~,idxL]                = min(abs(f-0.025));
[~,idxU]                = min(abs(f-0.975));

