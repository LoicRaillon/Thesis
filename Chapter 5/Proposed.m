% Proposed algorithm, section 5.3.4.2
% Adaptation of the recursive prediction error method into a sequential
% Monte Carlo framework
function out = Proposed(sys)

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
FIM                         = sys.FIM;                                      % Inverse approximation of Fisher information matrix
Np                          = sys.dim(1);                                   % Number of R,C parameters
Nx                          = sys.dim(2);                                   % Number of states
Nu                          = sys.dim(4);                                   % Number of inputs
Nprt                        = sys.Nprt;                                     % Number of particles
T                           = max(size(y));                                 % Number of samples
lb                          = sys.lb;                                       % Lower bounds
ub                          = sys.ub;                                       % Upper bounds
lbArray                     = bsxfun(@times,lb,ones(Np,Nprt));              % duplicate to array for easier indexing
ubArray                     = bsxfun(@times,ub,ones(Np,Nprt));              % ...
alpha                       = sys.alpha;                                    % percentage of parameter update

if ~isfield(sys,'ahold')                                                    % Check of sys.ahold exists
    ahold                   = 'zoh';                                        % if not, by default 'zoh'
else                                                                        % if specified
    ahold                   = sys.ahold;                                    % get user choice
end                                                                         % ...
if ~isfield(sys,'sig')                                                      % Check if sys.sig exists
    sig                     = [];                                           % if not, set sig to empty
    Q                       = sys.Q;                                        % get process noise covariance defined by user
else                                                                        % else
    sig                     = sys.sig;                                      % get incremental variance factor
end                                                                         % ...
R                           = sys.R;                                        % get measurement noise covariance

dx                          = zeros(Nx,Np,Nprt);                            % Allocation partial derivative states
dP                          = zeros(Nx,Nx,Np,Nprt);                         % Allocation partial derivative state covariance
xS                          = zeros(Nx+Np,T);                               % Allocation states and parameters posterior mean
G                           = zeros(Nx,2*Nu);                               % Allocation augemente input matrix
wlog                        = zeros(Nprt,1);                                % Allocation logarithm of importance weights
Ix                          = eye(Nx);                                      % Identity matrix [Nx x Nx]
Icte                        = eye(Np)*1e-10;                                % Ensure PSD inverse FIM
Onp                         = zeros(Nprt,Np);
afoh                        = zeros(Nu,T-1);                                % alpha (5.4), by default set to zero ('zoh')
if strcmp(ahold,'foh')                                                      % ...
    afoh                    = (u(:,2:end) - u(:,1:end-1))./dt;              % ...
end                                                                         % ...

xS(1:Nx,1)                  = mean(x,2);                                    % Save state posterior mean
xS(Nx+1:end,1)              = mean(m,2);                                    % Save parameter posterior mean
Tstr                        = num2str(T);
for k = 2:T
    disp([num2str(k), '/ ', Tstr])
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
            Q               = diag(diag(lyap(mdl.A,V)));                    % ...
        end                                                                 % ...
        S                   = abs(H*(F*P(:,:,i)*F' + Q)*H' + R);            % one step prediction error covariance (5.14)
        e                   = y(k)-H*(F*x(:,i)+G(:,1:Nu)*(afoh(:,k-1)*dt... % one step prediction error (5.15)
                              +u(:,k-1))-G(:,Nu+1:2*Nu)*afoh(:,k-1)); % ...                
        wlog(i)             = -0.5*(log(S) + e'./S.*e);                     % log of importance weights (5.30)
    end
    wp                      = exp(wlog - max(wlog));                        % Normalizing
    wp                      = wp/sum(wp);                                   % ...
    idx                     = sysresample(wp);                              % systematic resampling at each time step, the particles are uniformly distributed after resampling, w_k-1 not necessary in (5.30)
    [~,ridx,rpt]            = unique(idx);                                  % find indexing of resampled particles
    x                       = x(:,idx);                                     % Resample
    P                       = P(:,:,idx);                                   % ...
    m                       = m(:,idx);                                     % ...
    dx                      = dx(:,:,idx);                                  % ...
    dP                      = dP(:,:,:,idx);                                % ...
    FIM                     = FIM(:,:,idx);                                 % ...
    sys.y                   = y(k);                                         % Pass current data for gradient evaluation
    sys.u                   = u(:,k-1);                                     % ...
    mp                      = m;
    Ngrad                   = rpt(end);
    gam                     = ones(Np,Ngrad);
    for i = 1:Ngrad                                                         % for each distinct particles
        j                   = ridx(i);
        sys.x               = x(:,j);                                       % Pass to the function  
        sys.P               = P(:,:,j);                                     % ...
        sys.dx              = dx(:,:,j);                                    % ...
        sys.dP              = dP(:,:,:,j);                                  % ...
        sys.FIM             = FIM(:,:,j);                                   % ...
        sys.afoh            = afoh(:,k-1);                                  % ...
        out                 = GradProposed(m(:,j),sys);                     % Evaluate gradient
        x(:,i)              = out.x;                                        % Update
        P(:,:,i)            = out.P;                                        % ...
        dx(:,:,i)           = out.dx;                                       % ...
        dP(:,:,:,i)         = out.dP;                                       % ...
        FIM(:,:,i)          = out.FIM;                                      % ...
        di                  = out.dir;                                      % ...
        if isempty(alpha)
            m(:,i)          = m(:,j) - di;                                  % mean of (5.49), search direction
        else
            f               = m(:,j) - di;
            lbg             = m(:,j) - alpha.*m(:,j);
            ubg             = m(:,j) + alpha.*m(:,j);
            lf              = f < lbg;
            uf              = f > ubg;
            gam(lf,i)       = abs((lbg(lf)-m(lf,j))./di(lf));
            gam(uf,i)       = abs((ubg(uf)-m(uf,j))./di(uf));
            m(:,i)          = f;
            m(lf,i)         = lbg(lf);
            m(uf,i)         = ubg(uf);     
        end
    end
    m                       = m(:,rpt);                                     % duplicate
    x                       = x(:,rpt);                                     % ...
    P                       = P(:,:,rpt);                                   % ...    
    dx                      = dx(:,:,rpt);                                  % ...
    dP                      = dP(:,:,:,rpt);                                % ...
    FIM                     = FIM(:,:,rpt);                                 % ...
    gam                     = gam(:,rpt);
    xS(1:Nx,k)              = mean(x,2);                                    % Save posterior state mean
    xS(Nx+1:end,k)          = mean(m,2);                                    % Save posterior parameter mean
    m                       = m + gam.*mvnrnd(...
                                    Onp,mean(FIM(:,:,ridx),3)+Icte)';       % Propgate parameters (5.49)
    m(m < lb)               = lbArray(m < lb);                              % Clip parameters outside bounds
    m(m > ub)               = ubArray(m > ub);                              % ...
end
out.xS                      = xS;                                           % Output desired quantities
out.m                       = mp;                                           % ...
out.x                       = x;                                            % ...