% Monte Carlo simulations for evaluating the performances of the proposed
% algorithm vs the Liu and West particle filter
clear all;
close all;
clc
%% ------------------------------------------------------------------------
% Import Twin Houses experiment data available at 
% https://pure.strath.ac.uk/portal/en/datasets/twin-houses-empirical-validation-dataset-experiment-2(94559779-e781-4318-8842-80a2b1201668).html
% -------------------------------------------------------------------------
dt                          = 60*60;                                        % sampling time [s]
imp                         = importdata('Twin_house_exp2_60min.xlsx');     % Import data
cut                         = 465:840;                                      % ROLBS sequence
T                           = length(cut);                                  % Length of data
Tkitchen                    = imp.data(:,3);                                % Temperature kitchen [°C]
Tdoorway                    = imp.data(:,4);                                % Temperature doorway [°C]
Tparents                    = imp.data(:,5);                                % Temperature parent's bedroom [°C]
Tliving010                  = imp.data(:,6);                                % Temperature living room at 10 cm from the ground [°C]
Tliving110                  = imp.data(:,7);                                % Temperature living room at 110 cm from the ground [°C]
Tliving170                  = imp.data(:,10);                               % Temperature living room at 170 cm from the ground [°C]
Tliving                     = mean([Tliving010,Tliving110,Tliving170],2);   % Mean temperature of the living room  [°C]
Tcorridor                   = imp.data(:,11);                               % Temperature corridor [°C]
Tbath010                    = imp.data(:,12);                               % Temperature in the bathroom at 10 cm from the ground [°C]
Tbath110                    = imp.data(:,13);                               % Temperature in the bathroom at 110 cm from the ground [°C]
Tbath170                    = imp.data(:,14);                               % Temperature in the bathroom at 170 cm from the ground [°C]
Tbath                       = mean([Tbath010,Tbath110,Tbath170],2);         % Mean temperature of the bathroom [°C]
Tchild010                   = imp.data(:,15);                               % Temperature in the child room at 10 cm from the ground [°C]
Tchild110                   = imp.data(:,16);                               % Temperature in the child room at 110 cm from the ground [°C]
Tchild170                   = imp.data(:,17);                               % Temperature in the child room at 170 cm from the ground [°C]
Tchild                      = mean([Tchild010,Tchild110,Tchild170],2);      % Mean temperaure of the child room [°C]
Text                        = imp.data(:,54);                               % Outdoor air temperature [°C]
Qliving                     = imp.data(:,28);                               % Electric heater power in living room [W]
Qbath                       = imp.data(:,29);                               % Electric heater power in bathroom [W]
Qchild                      = imp.data(:,30);                               % Electric heater power in child room  [W]
Sr_east                     = imp.data(:,59);                               % solar radiation measured on the east vertical wall [W/m²]
Sr_south                    = imp.data(:,60);                               % solar radiation measured on the south vertical wall [W/m²]
Sr_west                     = imp.data(:,61);                               % solar radiation measured on the west vertical wall [W/m²]
%% ------------------------------------------------------------------------
% House geometry
% -------------------------------------------------------------------------
h                           = 2.60;                                         % heigth
Sliving_west                = 6.46*h;                                       % surface of the west wall of the living room 
Sbath_east                  = 3.31*h;                                       % surface of the east wall of the bathroom 
Schild_east                 = 2.88*h;                                       % surface of the east wall of the child room 
Sliving_south               = 5.20*h;                                       % surface of the south wall of the living room 
Schild_south                = 3.88*h;                                       % surface of the south wall of the child room
Swliving_west               = 1.29;                                         % surface of the west window of the living room 
Swliving_south              = 1.84+4.11;                                    % surface of the south windows of the living room
Swchild_south               = 1.29;                                         % surface of the south window of the child room
Swbath_east                 = 1.29;                                         % surface of the east window of the bethroom
Sgliving                    = 33.65;                                        % ground surface living room 
Sgchild                     = 11.19;                                        % ground surface child room
Sgbath                      = 6.92;                                         % ground surface bathroom
Sgcorridor                  = 5.46;                                         % ground surface corridor
Sgkitchen                   = 7.44;                                         % ground surface kitchen
Sgdoorway                   = 5.84;                                         % ground surface doorway
Sgparents                   = 11.19;                                        % ground surface parent's room
Sg1                         = Sgliving+ Sgchild + Sgbath + Sgcorridor;      % ground surface south zone
Sg2                         = Sgkitchen + Sgdoorway + Sgparents;            % ground surface north zone
Vspaces                     = [Sgliving;Sgchild;Sgbath;Sgcorridor]./Sg1;    % Normalized surfaces of the south zone
Vspaces2                    = [Sgkitchen;Sgdoorway;Sgparents]./Sg2;         % Normalized surfaces of the north zone
S1_wow_west                 = Sliving_west - Swliving_west;                 % south zone: west wall surface without windows
S1_wow_south                = Sliving_south + Schild_south ...              % south zone: south wall surface without windows
                                - Swliving_south - Swchild_south;           % ...
S1_wow_east                 = Schild_east + Sbath_east - Swbath_east;       % south zone: east wall surface without windows
%% ------------------------------------------------------------------------
% Inputs / Output
% -------------------------------------------------------------------------
Srw                         = S1_wow_west*Sr_west+S1_wow_east*Sr_east+...   % Total solar radiation on the walls [W]
                                                    S1_wow_south*Sr_south;  % ...
Sri                         = Swliving_west*Sr_west+...                     % Total solar radiation through the windows [W]
             (Swliving_south+Swchild_south)*Sr_south + Swbath_east*Sr_east; % ...
Qhvac1                      = Qliving + Qchild  + Qbath;                    % Electric power injected in the south zone [W]
Tzone1_all                  = [Tliving Tbath Tchild Tcorridor]*Vspaces;     % south zone temperature, weighted mean according to ground surfaces [°C]
Tzone2_all                  = [Tkitchen Tdoorway Tparents]*Vspaces2;        % north zone temperature, weighted mean according to ground surfaces [°C]
u                           = [Text(cut) Tzone2_all(cut) Srw(cut) ...       % inputs
                               Sri(cut) Qhvac1(cut)]';                      % ...
%% ------------------------------------------------------------------------
% Simulate Output
% -------------------------------------------------------------------------
sys.fun                     = @TwTaTm_aIRcoRoiRb;                           % model, Figure 5.1
sys.dt                      = dt;                                           % sampling time
sys.u                       = u;                                            % input
sys.ahold                   = 'foh';                                        % first order hold
sys.sc                      = 1e8;                                          % scaling of thermal capacities
sys.H                       = [0 1 0];                                      % output matrix
Np                          = 9;                                            % Number of R,C parameters
Nx                          = 3;                                            % Number of states
Ny                          = 1;                                            % Number of outputs
Nu                          = 5;                                            % Number of inputs
sys.dim                     = [Np Nx Ny Nu];                                % Pass dimensions
%                             [ lb      theta    ub ]
Ro                          = [1e-4 ; 4.85e-2 ; 2e-1];
Ri                          = [1e-4 ; 2.18e-3 ; 2e-1];
Rm                          = [1e-4 ; 2.67e-3 ; 2e-1];
Rb                          = [1e-4 ; 5.07e-3 ; 2e-1];
Cw                          = [1e-4 ; 3.38e-2 ;  1e1];
Ca                          = [1e-4 ; 7.47e-3 ;  1e1];
Cm                          = [1e-4 ; 1.43e-1 ;  1e1];
aI                          = [1e-4 ; 2.42e-1 ;  1e0];
Rso                         = [1e-4 ; 2.02e-3 ; 2e-1];
gather                      = [Ro Ri Rm Rb Cw Ca Cm aI Rso];                % gather in one variable
theta                       = gather(2,:)';                                 % get parameters
sys.lb                      = gather(1,:)';                                 % get lower bounds
sys.ub                      = gather(3,:)';                                 % get upper bounds
sys.alpha                   = [];                                           % 
sig                         = [];                                           % Don't use scaling of incremental variance
Q                           = diag([1e-4 1e-4 1e-4]);                       % process noise covariance
R                           = 1e-4;                                         % measurement noise covariance
sys.Q                       = Q;                                            % Pass to structure
sys.R                       = R;                                            % ...
sys.delta                   = 0.99;                                         % Discount factor Lui and West PF
mdl                         = feval(sys.fun,theta,sys);                     % Model evaluation
F                           = expm(mdl.A*dt);                               % Discretization (5.3)
bis                         = F - eye(Nx);                                  % ...
G                           = zeros(Nx,2*Nu);                               % ...
G(:,1:Nu)                   = mdl.A\bis*mdl.B;                              % ... Bd0
afoh                        = zeros(Nu,T-1);                                % alpha (5.4), by default set to zero ('zoh')
if strcmp(sys.ahold,'foh')                                                  % If first order hold
    G(:,Nu+1:2*Nu)          = mdl.A\(-mdl.A\bis + F*dt)*mdl.B;              % ... Bd1
    afoh                    = (u(:,2:end) - u(:,1:end-1))./dt;              % ...
end                                                                         % ...
if ~isempty(sig)                                                            % Process noise covariance
    Cs                      = sig*sig';                                     % solve Lyapunov equation
    V                       = Cs - F*Cs*F';                                 % ...
    Q                       = lyap(A,V);                                    % ...
    sys.sig                 = sig;                                          % Overwrite sys.Q if sig not empty
end                                                                         % ...
Qstd                        = chol(Q,'lower');                              % standard deviations
Rstd                        = chol(R,'lower');                              % ...
x0                          = [27 28 27]';                                  % initial states
P0                          = diag([0.5 0.1 0.5]);                          % initial states covariance
x0std                       = sqrt(diag(P0));                               % initial states std
xs                          = zeros(Nx,T);                                  % Allocation simulated states
xs(:,1)                     = x0;                                           % ...
for k = 2:T                                                                 % simulate
    xs(:,k)                 = F*xs(:,k-1) + G(:,1:Nu)*(afoh(:,k-1)*dt ...   % ...
                              + u(:,k-1)) - G(:,Nu+1:2*Nu)*afoh(:,k-1) ...  % ...
                              + Qstd*randn(Nx,1);                           % ...
end                                                                         % ...
y                           = sys.H*xs + bsxfun(@times,Rstd,randn(Ny,T));   % ...
sys.y                       = y;                                            % Pass to structure
%% ------------------------------------------------------------------------
% Start Monte Carlo runs
% -------------------------------------------------------------------------
MC                          = 1;                                            % Number of Monte Carlo runs
Nprt1                       = 4e3;                                          % Number of particles (RPEM)
Nprt2                       = 12e3;                                         % Number of particles (Liu and West PF)
xS1                         = zeros(Nx+Np,T,MC);                            % Allocation state and parameter posterior mean for each run (RPEM)
mS1                         = zeros(Np,Nprt1,MC);                           % Allocation parameter posteriors (RPEM)
xS2                         = zeros(Nx+Np,T,MC);                            % Allocation state and posterior mean for each run (Liu and West PF)
mS2                         = zeros(Np,Nprt2,MC);                           % Allocation parameter posteriors (RPEM)
m0std                       = theta.*0.1;                                   % initial std of parameters
FIM                         = diag([1e-6 1e-7 1e-7 1e-7 1e-6 1e-7 ...       % Initial inverse of FIM 
                                                        1e-5 1e-5 1e-7]);   % ...
for mc = 1:MC 
    m0std2                  = 2.*m0std;                                     % parameters random init between +/- 2*std and 3*std
    m0std3                  = 3.*m0std;                                     % ...
    rdm                     = rand(Np,1);                                   % ...
    lbr                     = theta + m0std2;                               % ...
    ubr                     = theta + m0std3;                               % ...
    lbr(rdm < 0.5)          = theta(rdm < 0.5) - m0std3(rdm < 0.5);         % ...
    ubr(rdm < 0.5)          = theta(rdm < 0.5) - m0std2(rdm < 0.5);         % ...
    m0rdm                   = lbr + (ubr - lbr).*rand(Np,1);                % ...
    x0std2                  = 2.*x0std;                                     % states random init between +/- 2*std and 3*std
    x0std3                  = 3.*x0std;                                     % ...
    rdm                     = rand(Nx,1);                                   % ...
    lbr                     = x0 + x0std2;                                  % ...
    ubr                     = x0 + x0std3;                                  % ...
    lbr(rdm < 0.5)          = x0(rdm < 0.5) - x0std3(rdm < 0.5);            % ...
    ubr(rdm < 0.5)          = x0(rdm < 0.5) - x0std2(rdm < 0.5);            % ...
    x0rdm                   = lbr + (ubr - lbr).*rand(Nx,1);                % ...
    % ------------------------------------------------------------------- %
    sys.Nprt                = Nprt1;                                        % Number of particles
    sys.FIM                 = repmat(FIM, [1 1 sys.Nprt]);                  % Re-Initialization
    sys.P                   = repmat(P0, [1 1 sys.Nprt]);                   % ...
    sys.x                   = bsxfun(@plus,x0rdm, ...                       % generate state particles
                                  bsxfun(@times,x0std,randn(Nx,sys.Nprt))); % ...
    sys.m                   = bsxfun(@plus,m0rdm, ...                       % generate parameter particles
                                  bsxfun(@times,m0std,randn(Np,sys.Nprt))); % ...
    out1                    = Proposed(sys);                                % Proposed algorithm
    xS1(:,:,mc)             = out1.xS;                                      % Save state and parameter posterior mean for each run
    mS1(:,:,mc)             = out1.m;                                       % Save parameter posteriors, particles 
    % ------------------------------------------------------------------- %
    sys.Nprt                = Nprt2;                                        % Number of particles
    sys.P                   = repmat(P0, [1 1 sys.Nprt]);                   % ...
    sys.x                   = bsxfun(@plus,x0rdm, ...                       % generate state particles
                                  bsxfun(@times,x0std,randn(Nx,sys.Nprt))); % ...
    sys.m                   = bsxfun(@plus,m0rdm, ...                       % generate parameter particles
                                  bsxfun(@times,m0std,randn(Np,sys.Nprt))); % ...
    out2                    = KernelDensity(sys);                           % Liu and West Particle filter
    xS2(:,:,mc)             = out2.xS;                                      % Save state and parameter posterior mean for each run
    mS2(:,:,mc)             = out2.m;                                       % Save parameter posteriors, particles 

    disp(num2str(mc))
end