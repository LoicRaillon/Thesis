function out = TwTaTm_aIRcoRoiRb(par,sys)

Np          = sys.dim(1); % Number of parameters
Nx          = sys.dim(2); % Number of states
Nu          = sys.dim(4); % Number of inputs
sc          = sys.sc;
A           = zeros(Nx,Nx);
B           = zeros(Nx,Nu);
dA1         = zeros(Nx,Nx,Np);
dB1         = zeros(Nx,Nu,Np);

Ro          = par(1);
Ri          = par(2);
Rm          = par(3);
Rb          = par(4);
Cw          = par(5);
Ca          = par(6); 
Cm          = par(7); 
aI          = par(8); 
Rco         = par(9);

A(1,1)      = -(Ro+Ri+Rco)/(Cw*Ri*(Ro+Rco)*sc);
A(1,2)      = 1/(Cw*Ri*sc);
A(2,1)      = 1/(Ca*Ri*sc);
A(2,2)      = -(Ri*Rm+Rb*Rm+Rb*Ri)/(Ca*Rb*Ri*Rm*sc);
A(2,3)      = 1/(Ca*Rm*sc);
A(3,2)      = 1/(Cm*Rm*sc);
A(3,3)      = -1/(Cm*Rm*sc);

B(1,1)      = 1/(Cw*(Ro+Rco)*sc);
B(1,3)      = Rco/(Cw*(Ro+Rco)*sc);
B(2,2)      = 1/(Ca*Rb*sc);
B(2,4)      = aI/(Ca*sc);
B(2,5)      = 1/(Ca*sc);

% Ro
dA1(1,1,1)  = 1/(Cw*(Ro+Rco)^2*sc);
dB1(1,1,1)  = -1/(Cw*(Ro+Rco)^2*sc);
dB1(1,3,1)  = -Rco/(Cw*(Ro+Rco)^2*sc);
% Ri
dA1(1,1,2)  = 1/(Cw*Ri^2*sc);
dA1(1,2,2)  = -1/(Cw*Ri^2*sc);
dA1(2,1,2)  = -1/(Ca*Ri^2*sc);
dA1(2,2,2)  = 1/(Ca*Ri^2*sc);
% Rm
dA1(2,2,3)  = 1/(Ca*Rm^2*sc);
dA1(2,3,3)  = -1/(Ca*Rm^2*sc);
dA1(3,2,3)  = -1/(Cm*Rm^2*sc);
dA1(3,3,3)  = 1/(Cm*Rm^2*sc);
% Rb
dA1(2,2,4)  = 1/(Ca*Rb^2*sc);
dB1(2,2,4)  = -1/(Ca*Rb^2*sc);
% Cw
dA1(1,1,5)  = (Ro+Ri+Rco)/(Cw^2*Ri*(Ro+Rco)*sc);
dA1(1,2,5)  = -1/(Cw^2*Ri*sc);
dB1(1,1,5)  = -1/(Cw^2*(Ro+Rco)*sc);
dB1(1,3,5)  = -Rco/(Cw^2*(Ro+Rco)*sc);
% Ca
dA1(2,1,6)  = -1/(Ca^2*Ri*sc);
dA1(2,2,6)  = (Ri*Rm+Rb*Rm+Rb*Ri)/(Ca^2*Rb*Ri*Rm*sc);
dA1(2,3,6)  = -1/(Ca^2*Rm*sc);
dB1(2,2,6)  = -1/(Ca^2*Rb*sc);
dB1(2,4,6)  = -aI/(Ca^2*sc);
dB1(2,5,6)  = -1/(Ca^2*sc);
% Cm
dA1(3,2,7)  = -1/(Cm^2*Rm*sc);
dA1(3,3,7)  = 1/(Cm^2*Rm*sc);
% aI
dB1(2,4,8)  = 1/(Ca*sc);
% Rco
dA1(1,1,9)  = 1/(Cw*(Ro+Rco)^2*sc);
dB1(1,1,9)  = -1/(Cw*(Ro+Rco)^2*sc);
dB1(1,3,9)  = Ro/(Cw*(Ro+Rco)^2*sc);

out.A       = A;
out.B       = B;
out.dA      = dA1;
out.dB      = dB1;
