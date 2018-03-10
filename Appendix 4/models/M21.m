function out = M21(par,sys)

Np              = sys.dim(1);
Nx              = sys.dim(2);
Nu              = sys.dim(4);

if isfield(sys,'fix')
    fix         = sys.fix;
    Nukwn       = sum(fix(1:Np));
else
    Nukwn       = 0;
end

Rw              = par(1);
Cw              = par(2);
Ci              = par(3);
aW              = par(4);

sc              = sys.sc;

A               = zeros(Nx,Nx);
A(1,1)          = -4/(Cw*Rw*sc);
A(1,2)          = 2/(Cw*Rw*sc);
A(2,1)          = 2/(Ci*Rw*sc);
A(2,2)          = -2/(Ci*Rw*sc);

B               = zeros(Nx,Nu);
B(1,1)          = 2/(Cw*Rw*sc);
B(1,2)          = aW/(Cw*sc);
B(2,3)          = 1/(Ci*sc);

out.dA          = [];
out.dB          = [];

if Nukwn ~= 0 
    
    dA1         = zeros(Nx,Nx,Np);
    dB1         = zeros(Nx,Nu,Np);

    % Rw
    dA1(1,1,1)  = 4/(Cw*Rw^2*sc);
    dA1(1,2,1)  = -2/(Cw*Rw^2*sc);
    dA1(2,1,1)  = -2/(Ci*Rw^2*sc);
    dA1(2,2,1)  = 2/(Ci*Rw^2*sc);
    dB1(1,1,1)  = -2/(Cw*Rw^2*sc);
    % Cw
    dA1(1,1,2)  = 4/(Cw^2*Rw*sc);
    dA1(1,2,2)  = -2/(Cw^2*Rw*sc);
    dB1(1,1,2)  = -2/(Cw^2*Rw*sc);
    dB1(1,2,2)  = -aW/(Cw^2*sc);
    % Ci
    dA1(2,1,3)  = -2/(Ci^2*Rw*sc);
    dA1(2,2,3)  = 2/(Ci^2*Rw*sc);
    dB1(2,3,3)  = -1/(Ci^2*sc);
    % aW
    dB1(1,2,4)  = 1/(Cw*sc);
    
    out.dA      = dA1(:,:,fix(1:Np));
    out.dB      = dB1(:,:,fix(1:Np));
end

out.A           = A;
out.B           = B;

