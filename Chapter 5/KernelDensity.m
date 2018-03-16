% Liu and West Particle filter, section 5.3.4.1
function out = KernelDensity(sys)

if ~isfield(sys,'Q') && ~isfield(sys,'sig')                                 % Safety
    disp('Scaling of the incremental variance or discrete process noise must be defined');
    return
end
y                           = sys.y;                                        % output y[k]
u                           = sys.u;                                        % inputs u[:,k]
H                           = sys.H;                                        % output matrix
dt                          = sys.dt;                                       % sampling time
x                           = sys.x;                                        % state particles
P                           = sys.P;                                        % state covariance
m                           = sys.m;                                        % parameter particles
Np                          = sys.dim(1);                                   % Number of R,C parameters
Nx                          = sys.dim(2);                                   % Number of states
Nu                          = sys.dim(4);                                   % Number of inputs
Nprt                        = sys.Nprt;                                     % Number of particles
T                           = max(size(y));                                 % Number of samples
lb                          = sys.lb;                                       % Lower bounds
ub                          = sys.ub;                                       % Upper bounds
lbArray                     = bsxfun(@times,lb,ones(Np,Nprt));              % duplicate to array for easier indexing
ubArray                     = bsxfun(@times,ub,ones(Np,Nprt));              % ...
delta                       = sys.delta;                                    % discount factor
h2                          = 1-((3*delta-1)/(2*delta))^2;                  % ...
a                           = sqrt(1-h2);                                   % ...
if ~isfield(sys,'ahold')                                                    % Check of sys.ahold exists
    ahold                   = 'zoh';                                        % if not, by default 'zoh'
else                                                                        % if specified
    ahold                   = sys.ahold;                                    % get user choice
end                                                                         % ...
afoh                        = zeros(Nu,T-1);                                % alpha (5.4), by default set to zero ('zoh')
if strcmp(ahold,'foh')                                                      % ...
    afoh                    = (u(:,2:end) - u(:,1:end-1))./dt;              % ...
end                                                                         % ...
if ~isfield(sys,'sig')                                                      % Check if sys.sig exists
    sig                     = [];                                           % if not, set sig to empty
    Q                       = sys.Q;                                        % get process noise covariance defined by user
else                                                                        % else
    sig                     = sys.sig;                                      % get incremental variance factor
end                                                                         % ...
R                           = sys.R;                                        % get measurement noise covariance

w                           = ones(Nprt,1)/Nprt;                            % Init importance weights
wlog                        = zeros(Nprt,1);                                % Allocation logarithm of importance weights
Ix                          = eye(Nx);                                      % Identity matrix [Nx x Nx]
xS                          = zeros(Nx+Np,T);                               % Allocation state and parameter posterior mean
xS(1:Nx,1)                  = mean(x,2);                                    % Save init state mean
xS(Nx+1:end,1)              = mean(m,2);                                    % Save init parameter mean
for k = 2:T
    xnm                     = m*w;                                          % Monte Carlo mean (5.35)
    tmp                     = bsxfun(@minus,m,xnm);                         % Monte Carlo variance (5.37)
    xnv                     = tmp*(bsxfun(@times,w,tmp'));                  % ...
    m                       = bsxfun(@plus,a*m,(1-a)*xnm);                  % kernel location (5.34)
    for i = 1:Nprt                                                          % For each particles
        mdl                 = feval(sys.fun,m(:,i),sys);                    % Model evaluation
        F                   = expm(mdl.A*dt);                               % Discretization (5.3)
        bis                 = F - Ix;                                       % ...
        G(:,1:Nu)           = mdl.A\bis*mdl.B;                              % ... Bd0
        if strcmp(ahold,'foh')                                              % If first order hold
            G(:,Nu+1:2*Nu)  = mdl.A\(-mdl.A\bis + F*dt)*mdl.B;              % ... Bd1
        end                                                                 % ...
        if ~isempty(sig)                                                    % Process noise covariance  
            Cs              = sig*sig';                                     % solve Lyapunov equation
            V               = Cs - F*Cs*F';                                 % ...
            Q               = lyap(A,V);                                    % ...
        end                                                                 % ...
        S                   = abs(H*(F*P(:,:,i)*F' + Q)*H' + R);            % one step prediction error covariance (5.14)
        e                   = y(k)-H*(F*x(:,i)+G(:,1:Nu)*(afoh(:,k-1)*dt... % one step prediction error (5.15)
                              +u(:,k-1))-G(:,Nu+1:2*Nu)*afoh(:,k-1));       % ...                
        wlog(i)             = -0.5*(log(S) + e'./S.*e);                     % log of importance weights (5.30)
    end
    wp                      = exp(wlog - max(wlog)).*w;                     % importance weights (5.30)    
    wp                      = wp/sum(wp);                                   % Normalizing (5.27)
    ess                     = 1/sum(wp.^2);                                 % effective sample size (5.28)
    if ess < 0.5*Nprt                                                       % if < to half the particles
        idx                 = sysresample(wp);                              % resample
        x                   = x(:,idx);                                     % ...
        P                   = P(:,:,idx);                                   % ...
        m                   = m(:,idx);                                     % ...
        w                   = ones(Nprt,1)/Nprt;                            % ...
    else                                                                    % else
        w                   = wp;                                           % save importance weights for next step
    end
    m                       = mvnrnd(m',h2.*xnv)';                          % Propgate parameters (5.33)
    m(m < lb)               = lbArray(m < lb);                              % Clip parameters outside bounds
    m(m > ub)               = ubArray(m > ub);                              % ...
    for i = 1:Nprt
        mdl                 = feval(sys.fun,m(:,i),sys);                    % Model evaluation
        F                   = expm(mdl.A*dt);                               % Discretization (5.3)
        bis                 = F - Ix;                                       % ...
        G(:,1:Nu)           = mdl.A\bis*mdl.B;                              % ... Bd0
        if strcmp(ahold,'foh')                                              % If first order hold
            G(:,Nu+1:2*Nu)  = mdl.A\(-mdl.A\bis + F*dt)*mdl.B;              % ... Bd1
        end                                                                 % ...
        if ~isempty(sig)                                                    % Process noise covariance  
            Cs              = sig*sig';                                     % solve Lyapunov equation
            V               = Cs - F*Cs*F';                                 % ...
            Q               = lyap(A,V);                                    % ...
        end                                                                 % ...
        x(:,i)              = F*x(:,i) + G(:,1:Nu)*(afoh(:,k-1)*dt ...      % Prior state mean (5.13)
                              + u(:,k-1)) - G(:,Nu+1:2*Nu)*afoh(:,k-1);     % ...                 
        P(:,:,i)            = F*P(:,:,i)*F' + Q;                            % Prior state covariance (5.13)
        S                   = H*P(:,:,i)*H' + R;                            % covariance one-step predictione error (5.14)
        K                   = P(:,:,i)*H'/S;                                % Kalman Gain (5.16)
        x(:,i)              = x(:,i) + K*(y(k) - H*x(:,i));                 % Posterior state mean (5.17)
        P(:,:,i)            = P(:,:,i)-K*S*K';                              % Posterior state covariance (5.17)
    end
    xS(1:Nx,k)              = x*w;                                          % Save state posterior mean
    xS(Nx+1:end,k)          = m*w;                                          % Save parameter posterior mean
end
out.xS                      = xS;                                           % output desired quantities
out.m                       = m;                                            % ...
out.x                       = x;                                            % ...