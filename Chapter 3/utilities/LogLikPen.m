% Compute penalysed LogLikelihood without gradient
% Only zero order hold, defined by default
function LL = LogLikPen(eta,sys)
y                               = sys.y;                                    % output 
u                               = sys.u;                                    % input
dt                              = sys.dt;                                   % sampling time
Np                              = sys.dim(1);                               % Number of R,C parameters
Nx                              = sys.dim(2);                               % Number of states
Ny                              = sys.dim(3);                               % Number of output
Nu                              = sys.dim(4);                               % Number of inputs
fix                             = sys.fix;                                  % logical indexing of fixed parameters
T                               = max(size(y));                             % Sample size
LL                              = 0.5*Ny*log(2*pi)*T;                       % Init loglikelihood
Oyx                             = zeros(Ny,Nx);                             % zeros [Ny x Nx]

etaB                            = eta(fix);                                 % unconstrained parameters
lbx                             = sys.lb(fix);                              % lower bounds
ubx                             = sys.ub(fix);                              % upper bounds
theta                           = eta;                                      % inverse transform to constrained parameters
[thetax,~,~,~]                  = inv_transform(etaB,'LowUp',lbx,ubx);      % ...
theta(fix)                      = thetax;                                   % ...
                                   
x                               = theta(Np+Nx+Ny+1:end);                    % initial  states
Psqrt                           = sys.Psqrt;                                % std of initial states
Qsqrt                           = diag(theta(Np+1:Np+Nx));                  % std of process noise
Rsqrt                           = diag(theta(Np+Nx+1:Np+Nx+Ny));            % std of measurement noise

model                           = feval(sys.fun,theta(1:Np),sys);           % Evaluate continuous model 
A                               = model.A;                                  % state matrix
B                               = model.B;                                  % input matrix
if isfield(model,'C')                                                       % if model.C exists
    C                           = model.C;                                  % use it
else                                                                        % else 
    C                           = sys.C;                                    % the output matrix must be specified in the structure
end
if isfield(model,'D')                                                       % if model.S exists    
    D                           = model.D;                                  % use it
else                                                                        % else
    D                           = zeros(Ny,Nu);                             % zero by default
end
F                               = expm(A*dt);                               % discretization
G                               = A\(F-eye(Nx))*B;                          % ...
%% ------------------------------------------------------------------------
% QR Kalman filter
% -------------------------------------------------------------------------
for k = 1:T
    Pre                         = [Rsqrt Oyx ; Psqrt*C' Psqrt];             % Pre array (4.16)
    [~,Post]                    = qr(Pre);                                  % post array (4.16)
    Ssqrt                       = Post(1:Ny,1:Ny);                          % std matrix residuals (4.16)
    nres                        = Ssqrt\(y(k)-C*x-D*u(:,k));                % standardized residuals (4.15)        
    nK                          = Post(1:Ny,Ny+1:Ny+Nx)';                   % Kalman gain, post array (4.16)
    Psqrt                       = Post(Ny+1:Ny+Nx,Ny+1:Ny+Nx);              % posterior std state matrix, post array (4.16)
    x                           = x + nK*nres;                              % posterior state mean (4.14)       
    LL                          = LL + 0.5*log(det(Ssqrt'*Ssqrt))...        % loglikelihood (4.28)
                                                 + 0.5*(nres')*nres;        % ...
    if k < T
        Pre                     = [Psqrt*F'; Qsqrt];                        % Pre array (4.23)
        [~,Post]                = qr(Pre);                                  % Post array (4.23)     
        Psqrt                   = Post(1:Nx,1:Nx);                          % prior std state matrix from post array (4.23)
        x                       = F*x + G*u(:,k);                           % Prior state mean (4.22)
    end
end

pen                             = 1e-4*sum(abs(lbx)./(thetax-lbx)+...       % penalty function 
                                                abs(ubx)./(ubx-thetax));    % ...
LL                              = LL + pen;                                 % Add penalty to negarive loglikelihood