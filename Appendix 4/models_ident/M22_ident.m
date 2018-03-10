function out = M22_ident(sys)

Np                  = 5;                                                    % Number of R,C parameters
Nx                  = 2;                                                    % Number of states
Ny                  = 1;                                                    % Number of outputs
Nu                  = 3;                                                    % Number of inputs
Ns                  = max(size(sys.y));                                     % Number of samples
idxu                = logical([1 0 1 1]);                                   % indexing inputs
sys.u               = sys.inputs(idxu,:);                                   % Get the inputs corresponding to the model
cte                 = sys.y(1);                                             % Fixed value
fix                 = [true(Np + Nx + Ny + (Nx-Ny),1) ; false];             % logical indexing of fixed parameter
Ng                  = sum(fix);                                             % Number of unknown parameters
sys.H               = [0 1];                                                % output matrix
sys.dim             = [Np Nx Ny Nu];                                        % vector of dimensions
sys.fix             = fix;                                                  % Pass indexing    
sys.fun             = @M22;                                                 % model function, M41
sys.Psqrt           = diag([0.5 0.1]);                                      % initial std of states

%                       lb     p     ub
Ro                  = [1e-5 ; 1e-2 ; 1];
Ri                  = [1e-5 ; 1e-2 ; 1];
Cw                  = [1e-5 ; 1e-1 ; 1e3];
Ci                  = [1e-5 ; 5e-3 ; 1e3];
aW                  = [1e-5 ; 1    ; 1e3];
sigw1               = [1e-8 ; 1e-2 ; 1];
sigw2               = [1e-8 ; 1e-2 ; 1];
sigv                = [1e-8 ; 1e-2 ; 1];
xw0                 = [10   ; 30   ; 50];
xi0                 = [10   ; cte  ; 50];
gather              = [Ro Ri Cw Ci aW sigw1 sigw2 sigv xw0 xi0];            % Gether in one variable
Vars                = {'Ro' 'Ri' 'Cw' 'Ci' 'aW' 'sigw1' ...                 % parameter names
                         'sigw2' 'sigv' 'xw0' 'xi0'}';                      % ...
lbx                 = gather(1,fix)';                                       % get lower bounds
sys.lb              = gather(1,:)';                                         % ...
ubx                 = gather(3,fix)';                                       % get upper bounds
sys.ub              = gather(3,:)';                                         % ...

eta                 = transform(gather(2,fix)','LowUp',lbx,ubx);
fun                 = @(x) SquareRootGradPen([x(1:9) ; cte], sys);
[phat,f,~,~,~,Hhat] = fminunc(fun,eta,sys.options);                         % Quasi Newton with BFGS
theta               = zeros(Np+2*Nx+Ny,1);                                  % Memory allocation for theta
[theta(fix),Jd]     = inv_transform(phat,'LowUp',lbx,ubx);                  % inverse transform to obtain theta
theta(~fix)         = cte;                                                  % output vector with fixed parameters
[~,dLL]             = SquareRootGrad(theta,sys);                            % get grad(eta)
g                   = Jd.*dLL;                                              % Chain rule gradient

J                   = diag((ubx-lbx)./((theta(fix)-lbx).*(ubx-theta(fix))));% Jacobian
Hess                = J*Hhat*J;                                             % Chain rule for Hessian(theta)
covM                = inv(Hess);                                            % Inverse Hessian
bis                 = diag(covM);                                           
if any(bis < 0)
    idx             = find(bis < 0);                                        % Check if var < 0                                      
    Vars{idx}                                                               % Display unidentifiable parameters
    disp('negative standard deviation !!! ')
    bis             = abs(bis);                                             % Prevent failure
end
se                  = sqrt(bis);                                            % ML standard deviations of parameters
p                   = 2*(1-tcdf(abs(theta(fix)./se),Ns - Ng));              % pvalue t-test
dpen                = Jd.*(1e-4*(abs(ubx)./(ubx-theta(fix)).^2-...          % Derivative of penalty function 
                            abs(lbx)./(theta(fix)-lbx).^2));                % ...

format shortE                                                               % 5 decimals
Table               = [Vars(fix) num2cell(theta(fix)) num2cell(se) ...      % Put results in Table
                        num2cell(p) num2cell(g) num2cell(dpen)];            % ....
T                   = cell2table(Table,'VariableNames',{'Parameters' ...    % Head of table
                                        'theta' 'se' 'p' 'grad' 'dpen'});   % ....
T                                                                           % display Table
    
out.theta           = theta;                                                % Thete MLE
out.LL              = -f;                                                   % Loglikelihood
out.grad            = g;                                                    % grad(eta)     
out.dpen            = dpen;                                                 % dpen(theta)
out.fix             = fix;                                                  % index of fixed parameters
out.Hhat            = Hess;                                                 % Hessian(theta)
out.Table           = T;                                                    % Table of results