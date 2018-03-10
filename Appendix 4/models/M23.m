function out = M23(par,sys)

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
Cw              = par(3);
Ci              = par(4);
aW              = par(5);
aI              = par(6);

sc              = sys.sc;

A               = zeros(Nx,Nx);
A(1,1)          = -(Ro+Ri)/(Cw*Ri*Ro*sc);
A(1,2)          = 1/(Cw*Ri*sc);
A(2,1)          = 1/(Ci*Ri*sc);
A(2,2)          = -1/(Ci*Ri*sc);

B               = zeros(Nx,Nu);
B(1,1)          = 1/(Cw*Ro*sc);
B(1,2)          = aW/(Cw*sc);
B(2,2)          = aI/(Ci*sc);
B(2,3)          = 1/(Ci*sc);

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
    % Cw
    dA1(1,1,3)  = (Ro+Ri)/(Cw^2*Ri*Ro*sc);
    dA1(1,2,3)  = -1/(Cw^2*Ri*sc);
    dB1(1,1,3)  = -1/(Cw^2*Ro*sc);
    dB1(1,2,3)  = -aW/(Cw^2*sc);
    % Ci
    dA1(2,1,4)  = -1/(Ci^2*Ri*sc);
    dA1(2,2,4)  = 1/(Ci^2*Ri*sc);
    dB1(2,2,4)  = -aI/(Ci^2*sc);
    dB1(2,3,4)  = -1/(Ci^2*sc);
    % aW
    dB1(1,2,5)  = 1/(Cw*sc);
    % aI
    dB1(2,2,6)  = 1/(Ci*sc);
    
    out.dA      = dA1(:,:,fix(1:Np));
    out.dB      = dB1(:,:,fix(1:Np));
end

out.A           = A;
out.B           = B;

