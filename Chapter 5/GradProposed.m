function out = GradProposed(theta,sys)

y                           = sys.y;                                        % output y[k]
u                           = sys.u;                                        % inputs u[:,k]
H                           = sys.H;                                        % output matrix
dt                          = sys.dt;                                       % sampling time
x                           = sys.x;                                        % states
dx                          = sys.dx;                                       % partial derivative states
P                           = sys.P;                                        % state covariance
dP                          = sys.dP;                                       % partial derivative state covariance
FIM                         = sys.FIM;                                      % Inverse approximation of Fisher information matrix
Np                          = sys.dim(1);                                   % Number of R,C parameters
Nx                          = sys.dim(2);                                   % Number of states
Ny                          = sys.dim(3);                                   % Number of outputs
Nu                          = sys.dim(4);                                   % Number of inputs
ahold                       = sys.ahold;                                    % discretization assumption: 'zoh' or 'foh'
afoh                        = sys.afoh;                                     % piecewise linear approximation between inputs
if ~isfield(sys,'sig')                                                      % Check if sys.sig exists
    sig                     = [];                                           % if not, set sig to empty
    Q                       = sys.Q;                                        % get process noise covariance defined by user
else                                                                        % else
    sig                     = sys.sig;                                      % get incremental variance factor
end                                                                         % ...
R                           = sys.R;                                        % Measurement noise covariance

Fmn                         = zeros(2*Nx,2*Nx,Np);                          % Allocation augmented state matrix, AIM (5.44) 
GG1                         = zeros(2*Nx,Nu,Np);                            % Allocation augmented input matrix, Bd0 (5.44)
GG2                         = zeros(2*Nx,Nu,Np);                            % Allocation augmented input matrix, Bd1 (5.44)
dS                          = zeros(Ny,Np);                                 % Allocation partial derivative one-step prediction error covariance
I2x                         = eye(2*Nx);                                    % Identity matrix [2*Nx x 2*Nx]
Oxx                         = zeros(Nx,Nx);                                 % Zero matrix [Nx x Nx]
idx1                        = false(2*Nx,1);                                % Logical indexing 
idx1(1:Nx)                  = true;                                         % ...
idx2                        = false(2*Nx,1);                                % Logical indexing 
idx2(Nx+1:end)              = true;                                         % ...

mdl                         = feval(sys.fun,theta,sys);                     % Model evaluation 
A                           = mdl.A;                                        % state matrix
B                           = mdl.B;                                        % input matrix
dA                          = mdl.dA;                                       % partial derivative state matrix
dB                          = mdl.dB;                                       % partial derivative input matrix
for j = 1:Np                                                                % Discretization
    FF                      = [A Oxx ; dA(:,:,j) A];                        % AM (5.44)
    Fmn(:,:,j)              = expm(FF*dt);                                  % AIM (4.25)
    bis                     = Fmn(:,:,j) - I2x;                             % ...
    GG1(:,:,j)              = FF\bis*[B;dB(:,:,j)];                         % [Bd0 ; dBd0] (5.44)
    if strcmp(ahold,'foh')                                                  % If first order hold
        GG2(:,:,j)          = FF\(-FF\bis+Fmn(:,:,j)*dt)*[B;dB(:,:,j)];     % [Bd1 ; dBd1] (5.44)
    end                                                                     % ...
end                                                                         % ...
F                           = Fmn(idx1,idx1,1);                             % get state matrix from AIM (5.44)
G(:,1:Nu)                   = GG1(idx1,:,1);                                % get Bd0 (5.44)
G(:,Nu+1:2*Nu)              = GG2(idx1,:,2);                                % get Bd1 (5.44) | G = [Bd0 Bd1]
if ~isempty(sig)                                                            
    Cs                      = sig*sig';                                     % solve Lyapunov equation
    V                       = Cs - F*Cs*F';                                 % ...
    Q                       = diag(diag(lyap(mdl.A,V)));                    % ...
end                                                                         % ...
for j = 1:Np                                                                % partial derivative Kalman filter propagation (5.43)
    dx(:,j)                 = Fmn(idx2,:,j)*[x;dx(:,j)]+GG1(idx2,:,j)*...   % ...
                                            (afoh*dt+u)-GG2(idx2,:,j)*afoh; % ...
    dF                      = Fmn(idx2,idx1,j);                             %    
    dP(:,:,j)               = F*dP(:,:,j)*F'+dF*P*F'+F*P*dF';               % + dQ/dtheta
end                                                                         % ...    
x                           = F*x + G(:,1:Nu)*(afoh*dt+u) - ...             % Kalman filter propagation (5.13)
                                        G(:,Nu+1:2*Nu)*afoh;                % ...
P                           = F*P*F' + Q;                                   % ...
S                           = R + H*P*H';                                   % One step prediction error covariance (5.14)
HtiS                        = H'/S;
K                           = P*HtiS;                                       % Kalman gain (5.16)
e                           = y - H*x;                                      % one step prediction error (5.15)
x                           = x + K*e;                                      % posterior state mean (5.17)
SKt                         = S*K';
P                           = P - K*SKt;                                    % posterior state covariance (5.17)
viS                         = e'/S;
de                          = -H*dx;                                        % partial derivative one step prediction error (5.42)
W                           = S + de*FIM*de';                               % recursion for inverse FIM update
L                           = FIM*de'/W;                                    % ...
for j = 1:Np
    dS(j)                   = H*dP(:,:,j)*H';                               % + dR(:,:,j); partial derivative one step prediction erro covariance (5.42)
    dK                      = dP(:,:,j)*HtiS-P*HtiS*dS(j)/S;                % partial derivative Kalman gain (5.45)
    dx(:,j)                 = dx(:,j) + dK*e - K*H*dx(:,j);                 % partial derivative posterior state mean (5.46)
    dKSKt                   = dK*SKt; 
    dP(:,:,j)               = dP(:,:,j) - dKSKt - K*dS(j)*K' - dKSKt';      % partial derivative posterior state covariance (5.46)
end
Ps                          = FIM-L*W*L';                                   % Recursion for parameter update
phi                         = S\dS;                                         % trace for non scalar cases
Ss                          = 2+phi*Ps*phi';                                % ...
Ls                          = (Ps*phi')/Ss;                                 % ...
FIM                         = Ps - Ls*Ss*Ls';                               % ...
dir                         = L*e + 0.5*FIM*(phi - viS*dS*viS')';           % search direction

out.x                       = x;                                            % output quantities required 
out.P                       = P;                                            % ...
out.dx                      = dx;                                           % ...
out.dP                      = dP;                                           % ...
out.dir                     = dir;                                          % ...
out.FIM                     = FIM;                                          % ...