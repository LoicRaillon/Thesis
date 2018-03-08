% Maximum likelihood estimation, Table 4.5 
% Quasi-Newton with penalty function and BFGS update

close all;
clc
addpath(userpath, 'models');
addpath(userpath, 'utilities');
addpath(userpath, 'results');
%% ------------------------------------------------------------------------
% Import Twin Houses experiment data available at 
% https://pure.strath.ac.uk/portal/en/datasets/twin-houses-empirical-validation-dataset-experiment-2(94559779-e781-4318-8842-80a2b1201668).html
% -------------------------------------------------------------------------
dt                          = 10*60;                                        % sampling time [s]
imp                         = importdata('Twin_house_exp2_10min.xlsx');     % Import data
imp2                        = load('M3_5alpha_chains.mat','p0S');           % Import initial parameters of MH
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
sys.Psqrt                   = diag([0.5 0.5 0.1]);                          % Initial standard deviations of states                         
options                     = optimoptions(@fminunc,...                     % optimization options fminunc
                                'Display','final-detailed',...              % ...
                                'Algorithm','quasi-newton',...              % ...
                                'GradObj','on',...                          % ...
                                'HessUpdate','bfgs',...                     % ...
                                'MaxIter',1e3,...                           % ...
                                'MaxFunEvals',1e4,...                       % ...
                                'TolX',1e-6);                               % ...
%                             [ lb  ;  theta0 ;  ub ]
Ro                          = [1e-3 ; 5.16e-2 ; 1e-1];
Ri                          = [1e-4 ; 2.20e-3 ; 5e-2];
Rm                          = [1e-4 ; 2.16e-3 ; 1e-1];
Rz                          = [1e-4 ; 5.06e-3 ; 1e-1];
Cw                          = [1e-3 ; 2.83e-2 ; 5e-1];
Cm                          = [1e-2 ; 1.44e-1 ; 5];
Ci                          = [1e-4 ; 6.68e-3 ; 1e-1];
awW                         = [1e-3 ; 1.80e-2 ; 1e-1];
awS                         = [1e-3 ; 8.66e-2 ; 1e-1];
aiW                         = [1e-2 ; 3.82e-1 ; 1];
aiE                         = [1e-2 ; 5.90e-1 ; 1];
aiS                         = [1e-2 ; 1.46e-1 ; 1];
sigw1                       = [1e-8 ; 8.30e-2 ; 1];
sigw2                       = [1e-8 ; 1.46e-2 ; 1];
sigw3                       = [1e-8 ; 1.15e-2 ; 1];
sigv                        = [1e-8 ; 1.79e-2 ; 1];
xw0                         = [20   ; 2.88e1  ; 40];
xm0                         = [20   ; 2.97e1  ; 40];
xi0                         = [20   ; y(1)    ; 40];
gather                      = [Ro Ri Rm Rz Cw Cm Ci awW awS aiW aiE aiS ... % gather in one variable
                                sigw1 sigw2 sigw3 sigv xw0 xm0 xi0];
sys.lb                      = gather(1,:)';                                 % get lower bounds
sys.ub                      = gather(3,:)';                                 % get upper bounds
eta                         = imp2.p0S(:,[1 2 3 5 6]);                      % Import initial parameter estimate MH
phatS                       = zeros(Nukwn,5);                               % Allocation parameter estimate
fS                          = zeros(1,5);                                   % Allocation loglikelihood
gS                          = zeros(Nukwn,5);                               % Allocation gradient
HhatS                       = zeros(Nukwn,Nukwn,5);                         % Allocation BFGS estimate
for i = 1:5
    fun                     = @(x) SquareRootGradPen([x(1:18) ; y(1)], sys);% Penalysed loglikelihood
    [phat,f,~,~,g,Hhat]     = fminunc(fun,eta(:,i),options);                % Quasi Newton with BFGS
    phatS(:,i)              = phat;
    fS(:,i)                 = f;
    gS(:,i)                 = g;
    HhatS(:,:,i)            = Hhat;
end

% save MLE_rdm_init
