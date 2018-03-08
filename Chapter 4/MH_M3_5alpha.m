close all;
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');
%% ------------------------------------------------------------------------
% Import Twin Houses experiment data available at 
% https://pure.strath.ac.uk/portal/en/datasets/twin-houses-empirical-validation-dataset-experiment-2(94559779-e781-4318-8842-80a2b1201668).html
% -------------------------------------------------------------------------
dt                          = 10*60;                                        % sampling time [s]
imp                         = importdata('Twin_house_exp2_10min.xlsx');     % Import data
cut                         = 2890:5000;                                    % ROLBS sequence
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
SrwWest                     = S1_wow_west*Sr_west;                          % South zone, solar radiation on the west wall [W] 
SrwEast                     = S1_wow_east*Sr_east;                          % South zone, solar radiation on the east wall [W]
SrwSouth                    = S1_wow_south*Sr_south;                        % South zone, solar radiation on the south wall [W]
SriWest                     = Swliving_west*Sr_west;                        % South zone, solar radiation through the west window [W]
SriEast                     = Swbath_east*Sr_east;                          % South zone, solar radiation through the east window [W]
SriSouth                    = (Swliving_south+Swchild_south)*Sr_south;      % South zone, solar radiation through the south windows [W]
Qhvac1                      = Qliving + Qchild  + Qbath;                    % Electric power injected in the south zone [W]
Tzone1_all                  = [Tliving Tbath Tchild Tcorridor]*Vspaces;     % south zone temperature, weighted mean according to ground surfaces [°C]
Tzone2_all                  = [Tkitchen Tdoorway Tparents]*Vspaces2;        % north zone temperature, weighted mean according to ground surfaces [°C]
y                           = Tzone1_all(cut);                              % output
u                           = [Text(cut) Tzone2_all(cut) Sr_west(cut) ...   % inputs
                                Sr_east(cut) Sr_south(cut) Qhvac1(cut)]';   % ...
                            
%% ------------------------------------------------------------------------
% Information on the problem / model 
% -------------------------------------------------------------------------
Np                          = 12;                                           % Number of R,C parameters
Nx                          = 3;                                            % Number of states
Ny                          = 1;                                            % Number of output
Nu                          = 6;                                            % Number of inputs 
% The parameters are organized as follow: [Np(RC) Nx(Process std) Ny(Measurement std) Nx(initial states)] 
fix                         = [true(Np + Nx + Ny + (Nx-Ny),1) ; false];     % Logical vector, false = fixed
Nukwn                       = sum(fix);                                     % number of unknown
sys.fun                     = @M3_5alpha;                                   % Model function
sys.H                       = [0 0 1];                                      % Output matrix
sys.dim                     = [Np Nx Ny Nu];                                % Vector of dimensions 
sys.fix                     = fix;                                          % indicatior of fixed parameters
sys.ahold                   = 'foh';                                        % 'foh': first order hold, 'zoh': zero order hold
sys.sc                      = 1e8;                                          % scaling of the thermal capacities
sys.dt                      = dt;                                           % sampling time
sys.y                       = y;                                            % output
sys.u                       = u;                                            % input 
sys.epsi                    = eye(Nukwn)*0.2;                               % step length (4.32)
sys.MC                      = 4500;                                         % Number of iterations 
sys.Psqrt                   = diag([0.5 0.5 0.1]);                          % Initial standard deviations of states
nMC                         = 5;                                            % number of runs
p0S                         = zeros(Nukwn,nMC);                             % Allocation for random parameters init (eta0)
LL_ALL                      = zeros(nMC,sys.MC);                            % Allocation for log posterior
p_ALL                       = zeros(nMC,Nukwn,sys.MC);                      % Allocation for trace of the Markov chains
acc_ALL                     = zeros(nMC,sys.MC);                            % Allocation for the accpetace ratio
lb                          = [1e-3 1e-4 1e-4 1e-4 1e-3 1e-2 1e-4 ...       % lower bounds of[Ro Ri Rm Rz Cw Cm Ci ...
    1e-3 1e-3 1e-2 1e-2 1e-2 -Inf -Inf -Inf 20 20 20]';                     % ... awW awS aiW aiE aiS sigw1 sigw2 sigw3 sigv xw0 xm0
ub                          = [1e-1 5e-2 1e-1 1e-1 5e-1 5 1e-1 ...          % upper bounds of[Ro Ri Rm Rz Cw Cm Ci ...
    1e-1 1e-1 1 1 1 Inf Inf Inf 40 40 40]';                                 % ... awW awS aiW aiE aiS sigw1 sigw2 sigw3 sigv xw0 xm0
sys.lb                      = lb;                                           % lower bounds
sys.ub                      = ub;                                           % uppper bounds
idx1                        = false(1,Nukwn);                               % Logical indexing RC parameters
idx1(1:Np)                  = true;                                         % ...
idx2                        = false(1,Nukwn);                               % Logical indexing std parameters
idx2(Np+1:Np+Nx+Ny)         = true;                                         % ...
idx3                        = false(1,Nukwn);                               % Logical indexing initial states
idx3(Np+Nx+Ny+1:Nukwn)      = true;                                         % ...
eta                         = zeros(Np+2*Nx+Ny,1);                          % Allocation unconstrained paramters
eta(~fix)                   = y(1);                                         % Initial states corresponding to the measured state: fixed

%% ------------------------------------------------------------------------
% Second-order Metropolis-Hastings with nMC random initial conditions
% -------------------------------------------------------------------------
for i = 1:nMC   
    % random draw from Beta(2,2,lb,ub) for the RC parameters (section 4.2.2.3)
    eta(idx1)               = transform(lb(idx1) + (ub(idx1)-lb(idx1)).*betarnd(2,2,[Np,1]),'LowUp',lb(idx1),ub(idx1));
    % random draw from Gamma(2,0.03) for the std parameters (section 4.2.2.3)
    eta(idx2)               = transform(gamrnd(2,0.03,[Nx+Ny,1]),'Log'); 
    % random draw from Beta(2,2,lb,ub) for the RC parameters (section 4.2.2.3)
    eta(idx3)               = transform(lb(idx3) + (ub(idx3)-lb(idx3)).*betarnd(2,2,[Nx-1,1]),'LowUp',lb(idx3),ub(idx3));     
    p0S(:,i)                = eta(fix);                                     % Save initial unconstrained parameters
    [LLS,pS,accS]           = MetropolisHastings(eta,sys);                  % Second-order Metropolis-Hastings
    LL_ALL(i,:)             = LLS;                                          % Save trace of log posterior
    p_ALL(i,:,:)            = pS;                                           % Save trace of parameters
    acc_ALL(i,:)            = accS;                                         % Save acceptace rate
end

% save M3_5alpha_chain3456
% eval(['!shutdown -s -f -t ' num2str(1800)])
