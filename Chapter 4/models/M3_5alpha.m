% Model M3: 3 states, FIgure 4.1

% Input: 
%   - par: Parameters of the thermal networks
%   - sys.dim: Dimension of the system 
%   - sys.fix: logical vector indicating the fixed parameters
%   - sys.sc: Scale of the thermal capacities
% Outputs: 
%   - out.A: state matrix
%   - out.B: input matrix
%   - out.dA: Jacobian state matrix
%   - out.dB: Jacobian input matrix

function out = M3_5alpha(par,sys)

Np              = sys.dim(1); % Number of parameters
Nx              = sys.dim(2); % Number of states
Nu              = sys.dim(4); % Number of inputs

if isfield(sys,'fix')
    fix         = sys.fix;
    Nukwn       = sum(fix(1:Np));
else
    Nukwn       = 0;
end

Ro              = par(1);
Ri              = par(2);
Rm              = par(3);
Rz              = par(4);
Cw              = par(5);
Cm              = par(6);
Ci              = par(7);
awW             = par(8);
awS             = par(9);
aiW             = par(10);
aiE             = par(11);
aiS             = par(12);

% Walls and windows surfaces, section 4.2.1
SwW             = 1.5506e+01;
SwS             = 1.6368e+01;
SiW             = 1.2900e+00;
SiE             = 1.2900e+00;
SiS             = 7.2400e+00;
sc              = sys.sc;

A               = zeros(Nx,Nx);
A(1,1)          = -(Ro+Ri)/(Cw*Ri*Ro*sc);
A(1,3)          = 1/(Cw*Ri*sc);
A(2,2)          = -1/(Cm*Rm*sc);
A(2,3)          = 1/(Cm*Rm*sc);
A(3,1)          = 1/(Ci*Ri*sc);
A(3,2)          = 1/(Ci*Rm*sc);
A(3,3)          = -((Rm+Ri)*Rz+Ri*Rm)/(Ci*Ri*Rm*Rz*sc);

B               = zeros(Nx,Nu);
B(1,1)          = 1/(Cw*Ro*sc);
B(1,3)          = (SwW*awW)/(Cw*sc);
B(1,5)          = (SwS*awS)/(Cw*sc);
B(3,2)          = 1/(Ci*Rz*sc);
B(3,3)          = (SiW*aiW)/(Ci*sc);
B(3,4)          = (SiE*aiE)/(Ci*sc);
B(3,5)          = (SiS*aiS)/(Ci*sc);
B(3,6)          = 1/(Ci*sc);

out.dA          = [];
out.dB          = [];

if Nukwn ~= 0 
    
    dA1         = zeros(Nx,Nx,Np);
    dB1         = zeros(Nx,Nu,Np);

    % Ro
    dA1(1,1,1)  = 1/(Cw*Ro^2*sc);
    dB1(1,1,1)  = -dA1(1,1,1);
    % Ri
    dA1(1,1,2)  = 1/(Cw*Ri^2*sc);
    dA1(1,3,2)  = -dA1(1,1,2);
    dA1(3,3,2)  = 1/(Ci*Ri^2*sc);
    dA1(3,1,2)  = -dA1(3,3,2);
    % Rm
    dA1(2,2,3)  = 1/(Cm*Rm^2*sc);
    dA1(2,3,3)  = -dA1(2,2,3);
    dA1(3,3,3)  = 1/(Ci*Rm^2*sc);
    dA1(3,2,3)  = -dA1(3,3,3);
    % Rz
    dA1(3,3,4)  = 1/(Ci*Rz^2*sc);
    dB1(3,2,4)  = -dA1(3,3,4);
    % Cw
    dA1(1,1,5)  = (Ro+Ri)/(Cw^2*Ri*Ro*sc);
    dA1(1,3,5)  = -1/(Cw^2*Ri*sc);
    dB1(1,1,5)  = -1/(Cw^2*Ro*sc);
    dB1(1,3,5)  = -(SwW*awW)/(Cw^2*sc);
    dB1(1,5,5)  = -(SwS*awS)/(Cw^2*sc);
    % Cm
    dA1(2,2,6)  = 1/(Cm^2*Rm*sc);
    dA1(2,3,6)  = -dA1(2,2,6);
    % Ci
    dA1(3,1,7)  = -1/(Ci^2*Ri*sc);
    dA1(3,2,7)  = -1/(Ci^2*Rm*sc);
    dA1(3,3,7)  = (Rm*Rz+Ri*Rz+Ri*Rm)/(Ci^2*Ri*Rm*Rz*sc);
    dB1(3,2,7)  = -1/(Ci^2*Rz*sc);
    dB1(3,6,7)  = -1/(Ci^2*sc);
    dB1(3,3,7)  = SiW*aiW*dB1(3,6,7);
    dB1(3,4,7)  = SiE*aiE*dB1(3,6,7);
    dB1(3,5,7)  = SiS*aiS*dB1(3,6,7);
    % awW
    dB1(1,3,8)  = SwW/(Cw*sc);
    % awS
    dB1(1,5,9)  = SwS/(Cw*sc);
    % aiW
    dB1(3,3,10) = SiW/(Ci*sc);
    % aiE
    dB1(3,4,11) = SiE/(Ci*sc);
    % aiS
    dB1(3,5,12) = SiS/(Ci*sc);
    
    out.dA      = dA1(:,:,fix(1:Np));
    out.dB      = dB1(:,:,fix(1:Np));
end

out.A           = A;
out.B           = B;

