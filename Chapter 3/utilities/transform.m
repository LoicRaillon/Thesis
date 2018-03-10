% Transformation from constrained space to unconstrained space 
% eta = f(eta), section 4.2.2.2
% Inputs: 
%   - theta: Constrained parameters
%   - type: 
%        'Log' --> Log transformation (4.33) 
%        'LowUp' --> Bound transformation (4.34)
%   - lb: Lower bounds if type = 'LowUp'
%   - ub: Upper bounds if type = 'LowUp'
% Outputs: 
%   - eta: Unconstrained parameters

function eta = transform(theta,type,lb,ub)

if nargin < 2
    eta         = theta;
else
    switch type
        case 'Log' % theta [0,Inf[
            eta = log(theta);
        case 'LowUp' % theta [lb,ub]
            z   = (theta-lb)./(ub-lb);
            eta = log(z./(1-z));
        otherwise
            disp('A correct transformation must be specified')
    end
end