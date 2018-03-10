% Compute Akaike's information criterion (AIC) and Bayesian information
% crtiterion (BIC)
% inputs: 
%   - LL: LogLikelihood
%   - Np: Number of unknown
%   - Ns: Sample size
function [AIC,BIC] = AIC_BIC(LL,Np,Ns)

AIC = 2*Np - 2*LL;
BIC = log(Ns)*Np -2*LL;
