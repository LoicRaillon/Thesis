% ML estimatio model M4, 
% Display table 4.7
function out = M4_ident(sys)

Np                  = 14;                                                   % Number of R,C parameters
Nx                  = 4;                                                    % Number of states
Ny                  = 1;                                                    % Number of outputs
Nu                  = 6;                                                    % Number of inputs
Ns                  = max(size(sys.y));                                     % Number of samples
sys.u               = sys.inputs2;                                          % Get the inputs corresponding to the model
cte                 = sys.y(1);                                             % Fixed value
fix                 = [true(Np + Nx + Ny + (Nx-Ny),1) ; false];             % logical indexing of fixed parameter
fix(4)              = false;                                                % Rs fixed
fix(9)              = false;                                                % Cs fixed
fix(18)             = false;                                                % sigw44 fixed
Ng                  = sum(fix);                                             % Number of unknown parameters
sys.H               = [0 0 0 1];                                            % output matrix
sys.dim             = [Np Nx Ny Nu];                                        % vector of dimensions
sys.fix             = fix;                                                  % Pass indexing
sys.fun             = @M4_5alpha;                                           % model function, M4
sys.Psqrt           = diag([0.5 0.5 0.5 0.1]);                              % initial std of states
%                     [ lb     p     ub]                                    
Ro                  = [1e-5 ; 5e-02 ; 1];
Ri                  = [1e-5 ; 2e-03   ; 1];
Rm                  = [1e-5 ; 2e-03   ; 1];
Rs                  = [1e-5 ; 2.52e-3 ; 1];
Rz                  = [1e-5 ; 5e-03   ; 1];
Cw                  = [1e-5 ; 2e-02   ; 1e3];
Ci                  = [1e-5 ; 5e-03   ; 1e3];
Cm                  = [1e-5 ; 1e-01   ; 1e3];
Cs                  = [1e-8 ; 2.52e-4 ; 1e3];
awW                 = [1e-8 ; 2e-02   ; 1];
awS                 = [1e-8 ; 7e-02   ; 1];
aiW                 = [1e-5 ; 4e-01   ; 1];
aiE                 = [1e-5 ; 6e-01   ; 1];
aiS                 = [1e-5 ; 1e-01   ; 1];
sigw1               = [1e-8 ; 8e-02   ; 1];
sigw2               = [1e-8 ; 2e-02   ; 1];
sigw3               = [1e-8 ; 2e-02   ; 1];
sigw4               = [1e-8 ; 1e-6    ; 1];
sigv                = [1e-8 ; 1e-02   ; 1];
xw0                 = [10   ; 29      ; 50];
xm0                 = [10   ; 28      ; 50];
xi0                 = [10   ; 29      ; 50];
xs0                 = [10   ; cte     ; 50];

gather              = [Ro Ri Rm Rs Rz Cw Ci Cm Cs awW awS aiW aiE aiS ...   % Gether in one variable
                        sigw1 sigw2 sigw3 sigw4 sigv xw0 xm0 xi0 xs0];      % ...
Vars                = {'Ro' 'Ri' 'Rm' 'Rs' 'Rz' 'Cw' 'Ci' 'Cm' 'Cs' ...     % parameter names
                       'awW' 'awS' 'aiW' 'aiE' 'aiS' 'sigw1' 'sigw2' ...    % ...
                       'sigw3' 'sigw4' 'sigv' 'xw0' 'xm0' 'xi0' 'xs0'}';    % ...
lbx                 = gather(1,fix)';                                       % get lower bounds
sys.lb              = gather(1,:)';                                         % ...
ubx                 = gather(3,fix)';                                       % get upper bounds
sys.ub              = gather(3,:)';                                         % ...

eta                 = transform(gather(2,fix)','LowUp',lbx,ubx);           % transform to unconstrained space
fun                 = @(x) SquareRootGradPen([x(1:3);Rs(2);x(4:7);Cs(2);... % Loglikelihood with penalty function 
                                       x(8:15);sigw4(2);x(16:19);cte],sys); % ...
[phat,f,~,~,~,Hhat] = fminunc(fun,eta,sys.options);                         % Quasi Newton with BFGS
theta               = zeros(Np+2*Nx+Ny,1);                                  % Memory allocation for theta
[theta(fix),Jd]     = inv_transform(phat,'LowUp',lbx,ubx);                 % inverse transform to obtain theta
theta(~fix)         = [Rs(2) ; Cs(2) ; sigw4(2) ; cte];                     % output vector with fixed parameters
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