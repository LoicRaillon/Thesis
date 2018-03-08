% Second-order Metropolis-Hastings
% inputs:
%   eta: unconstrained parameter: sample space
%   sys: structure with data, constants, options, etc
% output
%   LLSave: Log posterior(eta) for the MC iterations
%   pSave: eta for the MC iterations
%   accRate: acceptance rate of the algorithm 

%% ------------------------------------------------------------------------
% Algorithm 1 section 4.2.2.2
% -------------------------------------------------------------------------
function [LLS,pS,accS] = MetropolisHastings(eta,sys)

MC                          = sys.MC;                                       % Number of iterations          
fix                         = sys.fix;                                      % Logical vector of unknown parameters
epsi                        = sys.epsi;                                     % Step length matrix
Nukwn                       = sum(fix);                                     % Number of unknown parameters
pS                          = zeros(Nukwn,MC);                              % Allocation for trace eta
LLS                         = zeros(1,MC);                                  % Allocation for trace log posterior
accS                        = zeros(1,MC);                                  % Allocation for trace accpetance ratio
count                       = 0;                                            % init accpetance ratio
eta1                        = eta;                                          % Pass fixed parameters
out                         = SquareRootGradHess(eta,sys);                  % Evaluation at eta(0) (Line 2)
LP                          = out.LP;                                       % Log posterior(eta(0))
dLP                         = out.dLP;                                      % Gradient(eta(0))
d2LP                        = out.d2LP;                                     % Hessian(eta(0))
E                           = epsi/d2LP;                                    % Covariance proposal distribution                                  
sig                         = chol(E,'lower');                              % Std proposal distribution
f                           = eta(fix) + 0.5.*E*dLP ;                       % Mean of proposal distribution 
for mc = 1:MC
    eta1(fix)               = f + sig*randn(Nukwn,1);                       % Suggest a new candidate (eta(star)) (Line 4)
    out                     = SquareRootGradHess(eta1,sys);                 % Evaluate candidates (Line 5)
    LP1                     = out.LP;                                       % Log posterior(eta(star))
    dLP1                    = out.dLP;                                      % Gradient(eta(star))
    d2LP1                   = out.d2LP;                                     % Hessian(eta(star))
    E1                      = epsi/d2LP1;
    sig1                    = chol(E1,'lower');
    f1                      = eta1(fix) + 0.5.*E1*dLP1;             
    prop                    = LogNormPdf(eta1(fix),f,sig);                  % q(theta(star)|theta(i-1)) (4.9)       
    prop1                   = LogNormPdf(eta(fix),f1,sig1);                 % q(theta(i-1)|theta(star)) (4.9)
    accp                    = min(1,exp(LP1 - LP + prop1 - prop));          % Acceptance probability (4.9) 
    if rand < accp                                                          % Accept suggested parameter? (Line 8)
        accS(mc)            = 1;                                            % If eta(star) accepted, acceptance rate(i) = 1
        count               = count + 1;                                    % count number of accepted parameters
        eta                 = eta1;                                         % eta(i) = eta(star)
        LP                  = LP1;                                          % Log posterior(eta(i)) = Log posterior(eta(star))
        dLP                 = dLP1;                                         % Gradient(eta(i)) = Gradient(eta(star))
        d2LP                = d2LP1;                                        % Hessian(eta(i)) = Hessian(eta(star))
        E                   = epsi/d2LP;
        sig                 = chol(E,'lower');
        f                   = eta(fix) + 0.5.*E*dLP ;                       % Mean of proposal distribution (Line 4)
    end
    pS(:,mc)                = eta(fix);                                     % Save eta
    LLS(mc)                 = LP;                                           % Save Log posterior distribution     
    disp(['ratio: ',num2str(mc), '/', num2str(count), ' | Log Posterior ', num2str(LP)]);
end

end                                                                         % End 2nd order MH

%% ------------------------------------------------------------------------
% Log normal probability function N(x|mu,S²)
% -------------------------------------------------------------------------
function y = LogNormPdf(x,mu,S)

if size(x,1) > size(x,2)
    x                       = x';
end

if size(mu,1) > size(mu,2)
    mu                      = mu';
end

logdetS = 2*sum(log(diag(S)));
ncoeff  = length(x)*log(2*pi) + logdetS;
e       = x - mu;
y       = -0.5*(ncoeff + e/S*e');

end
