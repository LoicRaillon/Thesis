% Compute detailed model M7
function [A,B] = model_detailed(R,C)

v           = [1;-1];
A1          = blkdiag(v,v,v,v,v,v,v,v);
v2          = [0;1;-1;0];
A2          = blkdiag(v2,v2,v2,v2);
v3          = [0;0;0;1];
A3          = [v3;v3;v3;v3];
A           = [A1 A2 A3];
G           = diag(1./R);
Nn          = length(C);                                                    % Nodes
Nb          = length(R);                                                    % Branches
Nz          = length(C(C==0));
Ns          = Nn - Nz;                                                      % States
u           = false(1,Nn+Nb);
utrue       = [1 5 9 13 17 19 21 23 29];
u(utrue)    = true;

K           = -A'*G*A;
K11         = K(1:Nz,1:Nz);
K12         = K(1:Nz,Nz+1:end);
K21         = K(Nz+1:end,1:Nz);
K22         = K(Nz+1:end,Nz+1:end);
Kb          = A'*G;
Kb1         = Kb(1:Nz,:);
Kb2         = Kb(Nz+1:end,:);
Cc          = diag(C(Nz+1:end));

A           = Cc\(-K21/K11*K12+K22);
Bs          = Cc\[-K21/K11*Kb1+Kb2 -K21/K11 eye(Ns)];
Bs2         = Bs(:,u);
B           = [sum(Bs2(:,1:4),2) Bs2(:,5:end)];
