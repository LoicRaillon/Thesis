% Transformation from unconstrained space to constrained space 
% theta = f^-1(eta), section 4.2.2.2
% Inputs: 
%   - eta: Unconstrained parameters
%   - type: 
%        'Log' --> Log transformation (4.33) 
%        'LowUp' --> Bound transformation (4.34)
%   - lb: Lower bounds if type = 'LowUp'
%   - ub: Upper bounds if type = 'LowUp'
% Outputs: 
%   - theta: Constrained parameters
%   - J: Jacobian, (4.36)
%   - dLogJ: First derivative of log(J)
%   - d2LogJ: Second derivative of log(J)
function [theta,J,dLogJ,d2LogJ] = inv_transform(eta,type,lb,ub)

if nargin < 2
    theta               = eta;
else
    switch type
        case 'Log' % theta [0,Inf[
            theta       = exp(eta);
            if nargout > 1
               J        = theta;  
               n        = length(theta);
               dLogJ    = ones(n,1);
               d2LogJ   = zeros(n,1);
            end
        case 'LowUp' % theta [lb,ub]
            e           = exp(-eta);
            theta       = lb + (ub-lb).*(1./(1+e));
            if nargout > 1
               J        = (e.*(ub-lb))./(1+e).^2;
               eeta     = exp(eta);
               dLogJ    = (1-eeta)./(1+eeta);
               d2LogJ   = (-2.*eeta)./(1+eeta).^2;
            end
        otherwise
            disp('A correct transformation must be specified')
    end
end