% inputs:
%   - eta: unconstrained parameter space: sample space
%   - sys: structure with data, constants, options, etc
% outputs
%   - out: structure 
%       - out.LP: Log posterior distribution(eta)
%       - out.dLP: gradient(eta)
%       - out.d2LP: Hessian(eta)

% The function Diff_QR belongs to:
%   - Kulikova, M. V. & Tsyganova, J. V., 2016. 
%       A unified square-root approach for the score and Fisher information matrix computation in linear dynamic systems.
%       Mathematics and Computers in Simulation, 119, pp.128–141. Available at: http://dx.doi.org/10.1016/j.matcom.2015.07.007.
%   - Tsyganova, Y. V. & Kulikova, M.V., 2012. 
%       On efficient parametric identification methods for linear discrete stochastic systems.
%       Automation and Remote Control, 73(6), pp.962–975. Available at: http://link.springer.com/10.1134/S0005117912060033.

function out = SquareRootGradHess(eta,sys)
%% ------------------------------------------------------------------------
% Get information
% -------------------------------------------------------------------------
y                               = sys.y;                                    % Output vector
u                               = sys.u;                                    % Input vector
H                               = sys.H;                                    % Output matrix
dt                              = sys.dt;                                   % Sampling time
Np                              = sys.dim(1);                               % Number of parameter
Nx                              = sys.dim(2);                               % Number of state
Ny                              = sys.dim(3);                               % Number of output
Nu                              = sys.dim(4);                               % Number of input
Nall                            = Np + Nx + Ny + Nx;                        % total size of the parameter vector
T                               = max(size(y));                             % Number of measurements
lb                              = sys.lb;                                   % Lower bounds
ub                              = sys.ub;                                   % Upper bounds
Psqrt                           = sys.Psqrt;                                % Standard deviation matrix of initial states
fix                             = sys.fix;                                  % Logical vector of unknwon parameters
Nukwn                           = sum(fix);                                 % Number of unknown parameters
idx1                            = false(Nall,1);                            % Logical indexing RC parameters
idx1(1:Np)                      = fix(1:Np);                                % ...
idx2                            = false(Nall,1);                            % Logical indexing std parameters
idx2(Np+1:Np+Nx+Ny)             = fix(Np+1:Np+Nx+Ny);                       % ...
idx3                            = false(Nall,1);                            % Logical indexing initial state 
idx3(Np+Nx+Ny+1:Nall)           = fix(Np+Nx+Ny+1:Np+Nx+Ny+Nx);              % ...
if ~isfield(sys,'ahold')                                                    % Check if sys.ahold is specified
    ahold                       = 'zoh';                                    % By default --> 'zoh' : zero order hold
else                                                                        % ...
    ahold                       = sys.ahold;                                % 'zoh' or 'foh' (first order hold)
end
if ~isfield(sys,'dH')                                                       % Check if sys.dH is specified                                         
    dH                          = zeros(Ny,Nx,Nukwn);                       % By default --> Jacobian of output matrix is zero
else                                                                        % ...
    dH                          = sys.dH;                                   % get the Jacobian of the output matrix if given
end

%% ------------------------------------------------------------------------
% Transform to constrained space, theta = f^-1(eta) | section 4.2.2.2
% -------------------------------------------------------------------------
theta                           = eta;
[theta(idx1),J1,dJ1,d2J1]       = inv_transform(eta(idx1),'LowUp',lb(idx1),ub(idx1));
[theta(idx2),J2,dJ2,d2J2]       = inv_transform(eta(idx2),'Log'); 
[theta(idx3),J3,dJ3,d2J3]       = inv_transform(eta(idx3),'LowUp',lb(idx3),ub(idx3));
                                    

%% ------------------------------------------------------------------------
% Evaluate continuous model given in sys.fun
% -------------------------------------------------------------------------
model                           = feval(sys.fun,theta(1:Np),sys);
A                               = model.A;                                  % state matrix
B                               = model.B;                                  % input matrix
dA                              = model.dA;                                 % Jacobian state matrix
dB                              = model.dB;                                 % Jacobian input matrix

%% ------------------------------------------------------------------------
% Memory anLLocation and initialization
% -------------------------------------------------------------------------
dx                              = zeros(Nx,Nukwn);                          % state derivative
dPsqrt                          = zeros(Nx,Nx,Nukwn);                       % state std matrix derivative
Oyx                             = zeros(Ny,Nx);                             % zeros [Ny x Nx]
Oxx                             = zeros(Nx,Nx);                             % zeros [Nx x Nx]
dres                            = zeros(Nukwn,1);                           % standardized residual derivative
dSsqrt                          = zeros(Nukwn,1);                           % std residuals derivative
dLL                             = zeros(Nukwn,1);                           % gradient of loglikelihood
d2LL                            = zeros(Nukwn,Nukwn);                       % Hessian estimate loglikelihood
ind22                           = false(2*Nx,1);                            % Logical indexing for time integration
ind22(Nx+1:end)                 = true;                                     % ...
afoh                            = zeros(Nu,T-1);                            % alpha (4.8), by default set to zero ('zoh')
dQsqrt                          = zeros(Nx,Nx,Nukwn);                       % partial derivative of std process noise matrix
dRsqrt                          = zeros(Ny,Ny,Nukwn);                       % partial derivative of std measurement noise matrix
Fmn                             = zeros(2*Nx,2*Nx,Nukwn);                   % Augmented state matrix (4.25)
GG1                             = zeros(2*Nx,Nu,Nukwn);                     % Augmented input matrix, Bd0 (4.25)
GG2                             = zeros(2*Nx,Nu,Nukwn);                     % Augmented input matrix, Bd1 (4.25)
LL                              = -0.5*Ny*log(2*pi)*T;                      % Init loglikelihood
I2x                             = eye(2*Nx);                                % Identity matrix

%% ------------------------------------------------------------------------
% Find unknown process standard deviations (indexing)
% -------------------------------------------------------------------------
Qsqrt                           = diag(theta(Np+1:Np+Nx));
if any(fix(Np+1:Np+Nx))
    c                           = sum(idx1);
    i                           = 1;
    for n = 1:Nx
        if fix(Np+n)
            dQsqrt(n,n,c+i)     = 1;
            i                   = i + 1;
        end
    end
end

%% ------------------------------------------------------------------------
% Find unknown measurement standard deviations (indexing)  
% -------------------------------------------------------------------------
Rsqrt                           = diag(theta(Np+Nx+1:Np+Nx+Ny));
if any(fix(Np+Nx+1:Np+Nx+Ny))
    c                           = sum(fix(1:Np+Nx));
    i                           = 1;
    for n = 1:Ny
        if fix(Np+Nx+n)
            dRsqrt(n,n,c+i)     = 1;
            i                   = i + 1;
        end
    end
end

%% ------------------------------------------------------------------------
% Find unknown initial states (indexing)  
% -------------------------------------------------------------------------
x                               = theta(Np+Nx+Ny+1:end);
if sum(fix(Np+Nx+Ny+1:end)) ~= 0
    c                           = sum(fix(1:Np+Nx+Ny));
    i                           = 1;
    for n = 1:Nx
        if fix(Np+Nx+Ny+n)
            dx(n,c+i)           = 1;
            i                   = i + 1;
        end
    end
end

%% ------------------------------------------------------------------------
% Integrate state/input matrices with their jacobians
% -------------------------------------------------------------------------
c                               = sum(idx1);
if c > 0
    for j = 1:c
        FF                      = [A Oxx ; dA(:,:,j) A];                    % AM (4.25)
        Fmn(:,:,j)              = expm(FF*dt);                              % AIM (4.25)
        bis                     = Fmn(:,:,j) - I2x;
        GG1(:,:,j)              = FF\bis*[B;dB(:,:,j)];                     % [Bd0 ; dBd0] (4.25)
        if strcmp(ahold,'foh')                                              
            GG2(:,:,j)          = FF\(-FF\bis+Fmn(:,:,j)*dt)*[B;dB(:,:,j)]; % [Bd1 ; dBd1] (4.25)
            afoh                = (u(:,2:end) - u(:,1:end-1))./dt;          % alpha (4.8)
        end
    end
end
F                               = Fmn(1:Nx,1:Nx,1);                         % get state matrix from AIM (4.25)
G(:,1:Nu)                       = GG1(1:Nx,:,1);                            % get Bd0 (4.25)
G(:,Nu+1:2*Nu)                  = GG2(1:Nx,:,2);                            % get Bd1 (4.25) | G = [Bd0 Bd1]

%% ------------------------------------------------------------------------
% Kalman recursion and sensitivity equations
% -------------------------------------------------------------------------
for k = 1:T
    for j = 1:Nukwn
        Pre                     = [Rsqrt Oyx ; Psqrt*H' Psqrt];             % Pre array (4.16)
        dPre                    = [dRsqrt(:,:,j) Oyx ; ...                  % Pre array (4.18)
            dPsqrt(:,:,j)*H'+Psqrt*dH(:,:,j)' dPsqrt(:,:,j)];               % ...
        [Post,dPost]            = Diff_QR(Pre,Ny+Nx,dPre);                  % Post array (4.16) and (4.18)
        Ssqrt                   = Post(1:Ny,1:Ny);                          % std matrix residuals (4.16) 
        dSsqrt(j)               = dPost(1:Ny,1:Ny);                         % derivative std matrix residuals (4.17)
        nres                    = Ssqrt\(y(k)-H*x);                         % standardized residuals (4.15) 
        dres(j)                 = -Ssqrt\(dSsqrt(j)*nres + ...              % derivative standardized residuals (4.21)
            dH(:,:,j)*x + H*dx(:,j) );                                      % ...
        dPsqrt(:,:,j)           = dPost(Ny+1:Ny+Nx,Ny+1:Ny+Nx);             % derivative psoterior std state matrix (4.17) 
        nK                      = Post(1:Ny,Ny+1:Ny+Nx)';                   % Kalman gain, post array (4.16) 
        dnK                     = dPost(1:Ny,Ny+1:Ny+Nx)';                  % derivative kalman gain (4.17) 
        dx(:,j)                 = dx(:,j) + dnK*nres + nK*dres(j);          % derivative posterior state mean (4.20)
        dLL(j)                  = dLL(j) - (trace(Ssqrt\dSsqrt(j)) ...      % derivative loglikelihood (4.29)
            + (nres')*dres(j));                                             % ...
    end
    Psqrt                       = Post(Ny+1:Ny+Nx,Ny+1:Ny+Nx);              % posterior std state matrix, post array (4.16)
    x                           = x + nK*nres;                              % posterior state mean (4.14) 
    LL                          = LL - 0.5*(log(det(Ssqrt'*Ssqrt)) + ...    % loglikelihood (4.28)
        (nres')*nres);                                                      % ...
    d1                          = Ssqrt.\dSsqrt;                            
    d2LL                        = d2LL + d1*d1' + dres*dres';               % Hessian estimate (4.30) 
    if k < T
        for j = 1:Nukwn     
            Pre                 = [Psqrt*F'; Qsqrt];                        % Pre array (4.23)
            dF                  = Fmn(Nx+1:2*Nx,1:Nx,j);                    % get partial derivative of state matrix w.r.t theta(j)  
            dPre                = [dPsqrt(:,:,j)*F'+Psqrt*dF';...           % Pre array (4.27)
                                   dQsqrt(:,:,j)];                          % ...
            [Post,dPost]        = Diff_QR(Pre,Nx,dPre);                     % Post array (4.23) and (4.27)
            dx(:,j)             = dF*x + F*dx(:,j) + ...                    % derivative prior state mean (4.24)
                            GG1(ind22,:,j)*(afoh(:,k)*dt + u(:,k)) - ...    % ...
                            GG2(ind22,:,j)*afoh(:,k);                       % ...
            dPsqrt(:,:,j)       = dPost(1:Nx,1:Nx);                         % derivative prior std state matrix (4.26)
        end
        x                       = F*x + G(:,1:Nu)*(afoh(:,k)*dt + u(:,k))...% Prior state mean (4.22)
                            - G(:,Nu+1:2*Nu)*afoh(:,k);                     % ...
        Psqrt                   = Post(1:Nx,1:Nx);                          % prior std state matrix from post array (4.23) 
    end
end

%% ------------------------------------------------------------------------
% Evaluate prior distributions 
% -------------------------------------------------------------------------
[pB,dpB,d2pB]                   = LogBeta(theta(idx1),2,2,lb(idx1),ub(idx1));    
[pG,dpG,d2pG]                   = LogGamma(theta(idx2),2,0.03); 
[pN,dpN,d2pN]                   = LogBeta(theta(idx3),2,2,lb(idx3),ub(idx3));    
Prior                           = sum(pB) + sum(pG) + sum(pN);              
dPrior                          = [dpB ; dpG ; dpN];
d2Prior                         = diag([d2pB ; d2pG ; d2pN]);
LP                              = LL + Prior;                               % log(posterior) = log(likelihood) + log(prior)

J                               = [J1 ; J2 ; J3];                           % Jacobian adjustment (4.36)
Jd                              = diag(J);                                  % ...      
Jacobian                        = sum(log(abs(J)));                         % log(|det(J)|)
dJ                              = [dJ1 ; dJ2 ; dJ3];                        % first derivative log(det(J))
d2J                             = [d2J1 ; d2J2 ; d2J3];                     % second derivative log(det(J))
out.LP                          = LP + Jacobian;                            % (4.35) 
out.dLP                         = Jd'*(dLL + dPrior) + dJ;                  % (4.37)
out.d2LP                        = Jd'*(d2LL - d2Prior)*Jd + diag(d2J);      % (4.37)

end

%% ------------------------------------------------------------------------
% rescaled log(Beta(x|a,b)) betweeb lower lb and upper bound ub
% -------------------------------------------------------------------------
% inputs:
%   - x: value to evaluate
%   - a,b: shape paramters
%   - lb, ub: lower and upper bounds
% outputs: 
%   - p: log(Beta)
%   - dp: dlog(Beta)
%   - d2p: d2log(Beta)
function [p,dp,d2p]             = LogBeta(x, a, b, lb, ub)    
    p                           = (a-1).*log(x-lb) + (b-1).*log(ub-x) - betaln(a,b) - (a+b-1).*log(ub-lb);
    dp                          = -(a-1)./(lb - x) - (b-1)./(ub-x);
    d2p                         = -(a-1)./(lb-x).^2 - (b-1)./(ub-x).^2;
end

%% ------------------------------------------------------------------------
% log(Gamma(x|a,b))
% -------------------------------------------------------------------------
% inputs:
%   - x: value to evaluate
%   - a: shape paramters
%   - b: expected values
% outputs: 
%   - p: log(Gamma)
%   - dp: dlog(Gamma)
%   - d2p: d2log(Gamma)
function [p,dp,d2p]             = LogGamma(x,a,b)
    p                           = (a-1).*log(x) - (x./b) - a.*log(b) - gammaln(a);
    dp                          = (a-1)./x - 1./b;
    d2p                         = -(a-1)./x.^2;
end
