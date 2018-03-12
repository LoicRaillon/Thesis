% Saltelli, A. et al (2010), 'Variance bases sensitivity analysis of model output. Desing and estimator
% for the total sensitivity index', Computer Physics Communications.
% Elsevier B.V 181(2), pp 259-270. doi: 10.1016/j.cpc.2009.09.018
%
% Thesis section 2.3 Factor prioritization (Algorithm page 65)
% inputs:
%   - tfm: model transformation
%   - par: output parameter
%   - num: number of input parameters
%   - ts: sampling time
%   - N: Number of iterations
%   - T: Similarity transformation if required in @fct, line 29
% output:
%   - S: FIrst order sensitivity indices
%   - ST: Total sensitvity indices
function [S,ST] = var_sens(tfm,par,num,ts,N,T)

yAll            = zeros(2*N,1);                                             % Allocation output model                                                       
Vi              = zeros(N,num);                                             % Allocation first order variance
VT              = Vi;                                                       % Allocaiton total variance
lb              = 5e-2;                                                     % lower bounds
ub              = 1-lb;                                                     % upper bounds
k               = 1;                                                        % count
for i = 1:N
    % a) Generate independent sample matrices
    T51         = LPTAU51(100+i,2*num);                                     % Quasi random generator, better distributed without the 100 first values 
    A           = lb + (ub - lb).*T51(1:num);                               % Rescaling between lb and ub to avoid numerical instabilities ...
    B           = lb + (ub - lb).*T51(num+1:2*num);                         % ... in the model transformations
    % b) Radial sampling
    Ab          = repelem(A,num,1);
    idx         = diag(true(num,1));
    Ab(idx)     = B;
    % c) Model evaluation
    X           = [A;B;Ab];
    y           = fct(tfm,par,ts,X,T);    
    yA          = y(1);
    yB          = y(2);
    yAb         = y(3:end);  
    yAll(k:k+1) = [yA ; yB]; % vector for total variance
    Vi(i,:)     = yB.*(yAb-yA);                                             % Numerator (2.40)
    VT(i,:)     = (yA-yAb).^2;                                              % ...
    k           = k+2;
end 
Vtot            = var(yAll);                                                % total variance
S               = mean(Vi,1)/Vtot;                                          % (2.40)
ST              = mean(VT,1)/2/Vtot;                                        % (2.40)

