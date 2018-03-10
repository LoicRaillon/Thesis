function out = M41(par,sys)

Np              = sys.dim(1);
Nx              = sys.dim(2);
Nu              = sys.dim(4);

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
aW              = par(10);
aI              = par(11);

sc              = sys.sc;

A               = zeros(Nx,Nx);
A(1,1)          = -(Ro+Ri)/(Cw*Ri*Ro*sc);
A(1,2)          = 1/(Cw*Ri*sc);
A(2,1)          = 1/(Ci*Ri*sc);
A(2,2)          = -(((Rm+Ri)*Rs+Ri*Rm)*Rz+Ri*Rm*Rs)/(Ci*Ri*Rm*Rs*Rz*sc);
A(2,3)          = 1/(Ci*Rm*sc);
A(2,4)          = 1/(Ci*Rs*sc);
A(3,2)          = 1/(Cm*Rm*sc);
A(3,3)          = -1/(Cm*Rm*sc);
A(4,2)          = 1/(Cs*Rs*sc);
A(4,4)          = -1/(Cs*Rs*sc);

B               = zeros(Nx,Nu);
B(1,1)          = 1/(Cw*Ro*sc);
B(1,3)          = aW/(Cw*sc);
B(2,2)          = 1/(Ci*Rz*sc);
B(2,3)          = aI/(Ci*sc);
B(2,4)          = 1/(Ci*sc);

out.dA          = [];
out.dB          = [];

if Nukwn ~= 0 
    
    dA1         = zeros(Nx,Nx,Np);
    dB1         = zeros(Nx,Nu,Np);

    % Ro
    dA1(1,1,1)  = 1/(Cw*Ro^2*sc);
    dB1(1,1,1)  = -1/(Cw*Ro^2*sc);
    % Ri
    dA1(1,1,2)  = 1/(Cw*Ri^2*sc);
    dA1(1,2,2)  = -1/(Cw*Ri^2*sc);
    dA1(2,1,2)  = -1/(Ci*Ri^2*sc);
    dA1(2,2,2)  = 1/(Ci*Ri^2*sc);
    % Rm
    dA1(2,2,3)  = 1/(Ci*Rm^2*sc);
    dA1(2,3,3)  = -1/(Ci*Rm^2*sc);
    dA1(3,2,3)  = -1/(Cm*Rm^2*sc);
    dA1(3,3,3)  = 1/(Cm*Rm^2*sc);
    % Rs
    dA1(2,2,4)  = 1/(Ci*Rs^2*sc);
    dA1(2,4,4)  = -1/(Ci*Rs^2*sc);
    dA1(4,2,4)  = -1/(Cs*Rs^2*sc);
    dA1(4,4,4)  = 1/(Cs*Rs^2*sc);
    % Rz
    dA1(2,2,5)  = 1/(Ci*Rz^2*sc);
    dB1(2,2,5)  = -1/(Ci*Rz^2*sc);
    % Cw
    dA1(1,1,6)  = (Ro+Ri)/(Cw^2*Ri*Ro*sc);
    dA1(1,2,6)  = -1/(Cw^2*Ri*sc);
    dB1(1,1,6)  = -1/(Cw^2*Ro*sc);
    dB1(1,3,6)  = -aW/(Cw^2*sc);
    % Ci
    dA1(2,1,7)  = -1/(Ci^2*Ri*sc);
    dA1(2,2,7)  = (Rm*Rs*Rz+Ri*Rs*Rz+Ri*Rm*Rz+Ri*Rm*Rs)/(Ci^2*Ri*Rm*Rs*Rz*sc);
    dA1(2,3,7)  = -1/(Ci^2*Rm*sc);
    dA1(2,4,7)  = -1/(Ci^2*Rs*sc);
    dB1(2,2,7)  = -1/(Ci^2*Rz*sc);
    dB1(2,3,7)  = -aI/(Ci^2*sc);
    dB1(2,4,7)  = -1/(Ci^2*sc);
    % Cm
    dA1(3,2,8)  = -1/(Cm^2*Rm*sc);
    dA1(3,3,8)  = 1/(Cm^2*Rm*sc);
    % Cs
    dA1(4,2,9)  = -1/(Cs^2*Rs*sc);
    dA1(4,4,9)  = 1/(Cs^2*Rs*sc);
    % aW
    dB1(1,3,10) = 1/(Cw*sc);
    % aI
    dB1(2,3,11) = 1/(Ci*sc);
    
    out.dA      = dA1(:,:,fix(1:Np));
    out.dB      = dB1(:,:,fix(1:Np));
end

out.A           = A;
out.B           = B;

