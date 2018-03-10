function out = M1(par,sys)

sc              = sys.sc;
Rso             = par(1);
Rwo             = par(2);
Rwi             = par(3);
Rsi             = par(4);
Cw              = par(5);
S               = 10;

A               = -(Rso+Rsi+Rwo+Rwi)/(Cw*(Rsi+Rwi)*(Rso+Rwo)*sc);
B               = [1/(Cw*(Rso+Rwo)*sc) 1/(Cw*(Rsi+Rwi)*sc)];
C               = -1/((Rsi+Rwi)*S);
D               = [0 1/((Rsi+Rwi)*S)];

out.A           = A;
out.B           = B;
out.C           = C;
out.D           = D;