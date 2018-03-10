% QR Kalman filter, save state estimate, loglikelihood and standardized residuals
function [e,xS,LL] = Residuals(theta,sys)
%% ------------------------------------------------------------------------
% Get information
% -------------------------------------------------------------------------
y                               = sys.y;                                    % Output vector
u                               = sys.u;                                    % Input vector
H                               = sys.H;                                    % Output matrix
dt                              = sys.dt;                                   % Sampling time
Np                              = sys.dim(1);                               % Number of parameters
Nx                              = sys.dim(2);                               % Number of states
Ny                              = sys.dim(3);                               % Number of outputs         
Nu                              = sys.dim(4);                               % Number of inputs   
T                               = max(size(y));                             % Length of data
x                               = theta(Np+Nx+Ny+1:end);                    % Initial state vector
Psqrt                           = sys.Psqrt;                                % Initial std state matrix
Qsqrt                           = diag(theta(Np+1:Np+Nx));                  % std process noise
Rsqrt                           = diag(theta(Np+Nx+1:Np+Nx+Ny));            % std measurement noise
if ~isfield(sys,'ahold')                                                    % Check if sys.ahold is specified
    ahold                       = 'zoh';                                    % By default --> 'zoh' : zero order hold
else                                                                        % ...
    ahold                       = sys.ahold;                                % 'zoh' or 'foh' (first order hold)
end

%% ------------------------------------------------------------------------
% Evaluate continuous model given in sys.fun
% -------------------------------------------------------------------------
model                           = feval(sys.fun,theta(1:Np),sys);
A                               = model.A;                                  % state matrix
B                               = model.B;                                  % Input matrix

%% ------------------------------------------------------------------------
% Memory anLLocation and initialization
% -------------------------------------------------------------------------
Oyx                             = zeros(Ny,Nx);                             % zeros [Ny x Nx]
afoh                            = zeros(Nu,T-1);                            % alpha (4.8), by default set to zero ('zoh')
e                               = zeros(T,Ny);                              % Allocation for standardized residuals
xS                              = zeros(Nx,T);                              % Allocation for state estimates
LL                              = zeros(Ny,T);                              % Allocation for log likelihood

%% ------------------------------------------------------------------------
% Discretization (4.7) and (4.8)
% -------------------------------------------------------------------------
F                               = expm(A*dt);
bis                             = F-eye(Nx);
G                               = zeros(Nx,2*Nu);
G(:,1:Nu)                       = A\bis*B;
if strcmp(ahold,'foh')
    G(:,Nu+1:2*Nu)              = A\(-A\bis + F*dt)*B;
    afoh                        = (u(:,2:end) - u(:,1:end-1))./dt;
end   

%% ------------------------------------------------------------------------
% QR Kalman filter
% -------------------------------------------------------------------------
for k = 1:T
    xS(:,k)                     = x;                                        % Save state estimate
    Pre                         = [Rsqrt Oyx ; Psqrt*H' Psqrt];             % Pre array (4.16)
    [~,Post]                    = qr(Pre);                                  % post array (4.16)
    Ssqrt                       = Post(1:Ny,1:Ny);                          % std matrix residuals (4.16)
    nres                        = Ssqrt\(y(k)-H*x);                         % standardized residuals (4.15)        
    nK                          = Post(1:Ny,Ny+1:Ny+Nx)';                   % Kalman gain, post array (4.16)
    Psqrt                       = Post(Ny+1:Ny+Nx,Ny+1:Ny+Nx);              % posterior std state matrix, post array (4.16)
    x                           = x + nK*nres;                              % posterior state mean (4.14)
    LL(k)                       = -0.5*Ny*log(2*pi) - ...                   % loglikelihood (4.28)
        0.5*log(det(Ssqrt'*Ssqrt)) - 0.5*(nres')*nres;                      % ...
    if k < T
        Pre                     = [Psqrt*F'; Qsqrt];                        % Pre array (4.23)
        [~,Post]                = qr(Pre);                                  % Post array (4.23)
        Psqrt                   = Post(1:Nx,1:Nx);                          % prior std state matrix from post array (4.23)
        x                       = F*x + G(:,1:Nu)*(afoh(:,k)*dt + u(:,k))...% Prior state mean (4.22)
            - G(:,Nu+1:2*Nu)*afoh(:,k);                                     % ...
    end
    e(k)                        = nres;                                     % Save standardized residuals
end
