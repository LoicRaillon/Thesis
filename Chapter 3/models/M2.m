function out = M2(par,sys)

sc              = sys.sc;

Ro              = par(1);
Ri              = par(2);
Cw              = par(3);
S               = 10;

A               = -(Ro+Ri)/(sc*Cw*Ri*Ro);
B               = [1/(Ro*sc*Cw) 1/(Ri*sc*Cw)];
C               = -1/(Ri*S);
D               = [0 1/(Ri*S)];

out.A           = A;
out.B           = B;
out.C           = C;
out.D           = D;