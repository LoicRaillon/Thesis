% Model M4: 4 states, FIgure 4.1

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

function out = M4_5alpha(par,sys)

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
Rs              = par(4);
Rz              = par(5);
Cw              = par(6);
Ci              = par(7);
Cm              = par(8);
Cs              = par(9);
awW             = par(10);
awS             = par(11);
aiW             = par(12);
aiE             = par(13);
aiS             = par(14);

% Walls and windows surfaces, section 4.2.1
SwW             = 1.5506e+01;
SwS             = 1.6368e+01;
SiW             = 1.2900e+00;
SiE             = 1.2900e+00;
SiS             = 7.2400e+00;
sc              = sys.sc;

A               = zeros(Nx,Nx);
A(1,1)          = -(Ro+Ri)/(Cw*Ri*Ro*sc);
A(1,2)          = 1/(Cw*Ri*sc);
A(2,1)          = 1/(Ci*Ri*sc);
A(2,2)          = -(((Rm+Ri)*Rs+Ri*Rm)*Rz+Ri*Rm*Rs)/(Ci*Ri*Rm*Rs*Rz*sc);
A(2,3)          = 1/(Ci*Rm*sc);
A(2,4)          = 1/(Ci*Rs*sc);
A(3,2)          = 1/(Cm*Rm*sc);
A(3,3)          = -A(3,2);
A(4,2)          = 1/(Cs*Rs*sc);
A(4,4)          = -A(4,2);

B               = zeros(Nx,Nu);
B(1,1)          = 1/(Cw*Ro*sc);
B(1,3)          = (SwW*awW)/(Cw*sc);
B(1,5)          = (SwS*awS)/(Cw*sc);
B(2,2)          = 1/(Ci*Rz*sc);
B(2,6)          = 1/(Ci*sc);
B(2,3)          = SiW*aiW*B(2,6);
B(2,4)          = SiE*aiE*B(2,6);
B(2,5)          = SiS*aiS*B(2,6);

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
    dA1(1,2,2)  = -dA1(1,1,2);
    dA1(2,2,2)  = 1/(Ci*Ri^2*sc);
    dA1(2,1,2)  = -dA1(2,2,2);
    % Rm
    dA1(2,2,3)  = 1/(Ci*Rm^2*sc);
    dA1(2,3,3)  = -dA1(2,2,3);
    dA1(3,3,3)  = 1/(Cm*Rm^2*sc);
    dA1(3,2,3)  = -dA1(3,3,3);
    % Rs
    dA1(2,2,4)  = 1/(Ci*Rs^2*sc);
    dA1(2,4,4)  = -dA1(2,2,4);
    dA1(4,4,4)  = 1/(Cs*Rs^2*sc);
    dA1(4,2,4)  = -dA1(4,4,4);
    % Rz
    dA1(2,2,5)  = 1/(Ci*Rz^2*sc);
    dB1(2,2,5)  = -dA1(2,2,5);
    % Cw
    dA1(1,1,6)  = (Ro+Ri)/(Cw^2*Ri*Ro*sc);
    dA1(1,2,6)  = -1/(Cw^2*Ri*sc);
    dB1(1,1,6)  = -1/(Cw^2*Ro*sc);
    dB1(1,3,6)  = -(SwW*awW)/(Cw^2*sc);
    dB1(1,5,6)  = -(SwS*awS)/(Cw^2*sc);
    % Ci
    dA1(2,1,7)  = -1/(Ci^2*Ri*sc);
    dA1(2,2,7)  = (Rm*Rs*Rz+Ri*Rs*Rz+Ri*Rm*Rz+Ri*Rm*Rs)/(Ci^2*Ri*Rm*Rs*Rz*sc);
    dA1(2,3,7)  = -1/(Ci^2*Rm*sc);
    dA1(2,4,7)  = -1/(Ci^2*Rs*sc);
    dB1(2,2,7)  = -1/(Ci^2*Rz*sc);
    dB1(2,3,7)  = -(SiW*aiW)/(Ci^2*sc);
    dB1(2,4,7)  = -(SiE*aiE)/(Ci^2*sc);
    dB1(2,5,7)  = -(SiS*aiS)/(Ci^2*sc);
    dB1(2,6,7)  = -1/(Ci^2*sc);
    % Cm
    dA1(3,3,8)  = 1/(Cm^2*Rm*sc);
    dA1(3,2,8)  = -dA1(3,3,8);
    % Cs
    dA1(4,4,9)  = 1/(Cs^2*Rs*sc);
    dA1(4,2,9)  = -dA1(4,4,9);
    % awW
    dB1(1,3,10) = SwW/(Cw*sc);
    % awS
    dB1(1,5,11) = SwS/(Cw*sc);
    % aiW
    dB1(2,3,12) = SiW/(Ci*sc);
    % aiE
    dB1(2,4,13) = SiE/(Ci*sc);
    % aiS
    dB1(2,5,14) = SiS/(Ci*sc);
    
    out.dA      = dA1(:,:,fix(1:Np));
    out.dB      = dB1(:,:,fix(1:Np));
end

out.A           = A;
out.B           = B;

