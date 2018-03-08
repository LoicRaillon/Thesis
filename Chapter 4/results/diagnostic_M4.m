% Compute:
%   - acceptance ratio
%   - Potential Scale Reduction Factor, R
%   - Effective sample size, neff
%   - Integrated autocorrelation time, IACT
%   - Time cosntants
%   - 95% probability mass of posterior distribution 

% Functions used: 
% - M4_5alpha: Model M4
% - inv_transform: theta = f^-1(eta) 
% - IACT_plot: Compute IACT and plot autocorrelation of the Markov Chains
% - psrf: Gelman-Rubin Potential Scale Reduction Factor, fonction accessible in the MCMC Diagnostics for Matlab at: http://becs.aalto.fi/en/research/bayes/mcmcdiag/  

clear all
close all
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');
addpath(userpath, 'export_fig');

%% -------------------------------------------------------------------------
% Import data for model M3
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

start                       = 501; % burnin
stop                        = sys.MC;
chain1                      = reshape(impM4_1.p_ALL(1,:,:),[Nukwn,stop]);
chain2                      = reshape(impM4_1.p_ALL(2,:,:),[Nukwn,stop]);
chain3                      = reshape(impM4_2.p_ALL(1,:,:),[Nukwn,stop]);
chain4                      = reshape(impM4_2.p_ALL(2,:,:),[Nukwn,stop]);
chain5                      = reshape(impM4_2.p_ALL(3,:,:),[Nukwn,stop]);

%% -------------------------------------------------------------------------
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
% Average acceptance ratio of the 5 Markov chains
% -------------------------------------------------------------------------
ratio1                      = sum(impM4_1.acc_ALL(:,start:stop),2)./length(impM4_1.acc_ALL(:,start:stop));
ratio2                      = sum(impM4_2.acc_ALL(:,start:stop),2)./length(impM4_2.acc_ALL(:,start:stop));
ratio                       = mean([ratio1 ; ratio2]);

%% ------------------------------------------------------------------------
% Convergence diagnosis, section 4.2.2.5
% -------------------------------------------------------------------------
% R : Potential scale factor reduction
% neff : Effective sample size

% First split each chains in two to have 2M chains
% Potential scale factor reduction
% Remove burn-in
% M chains
% N length chains
% Np number of parameters
% X = [NxNpxM] matrix 
% N = [MC - burnin]/2
N                           = (sys.MC - start+1)/2;
ch1                         = theta1_M4(:,start+1:start+N)';
ch2                         = theta1_M4(:,start+N:end)';
ch3                         = theta2_M4(:,start+1:start+N)';
ch4                         = theta2_M4(:,start+N:end)';
ch5                         = theta3_M4(:,start+1:start+N)';
ch6                         = theta3_M4(:,start+N:end)';
ch7                         = theta4_M4(:,start+1:start+N)';
ch8                         = theta4_M4(:,start+N:end)';
ch9                         = theta5_M4(:,start+1:start+N)';
ch10                        = theta5_M4(:,start+N:end)';
X                           = cat(3,ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9,ch10);
[R,neff,~,~,~,~,~]          = psrf(X);

%% ------------------------------------------------------------------------
% Integrated autocorrelation time (IACT), section 4.2.2.4
% -------------------------------------------------------------------------
nlags                       = 200; % number of lags
ncol                        = 3; % number of columns in the figure
iact1                       = IACT_plot(theta1_M4,start,nlags,ncol);
iact2                       = IACT_plot(theta2_M4,start,nlags,ncol);
iact3                       = IACT_plot(theta3_M4,start,nlags,ncol);
iact4                       = IACT_plot(theta4_M4,start,nlags,ncol);
iact5                       = IACT_plot(theta5_M4,start,nlags,ncol);
iact                        = [iact1 iact2 iact3 iact4 iact5];
mmiact                      = [min(iact,[],2) max(iact,[],2)]; % Min/Max IACT

