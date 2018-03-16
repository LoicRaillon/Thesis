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
y                           = Tzone1_all(cut);                              % output
%% ------------------------------------------------------------------------
% Sequential identification Twin Houses Experiment
% -------------------------------------------------------------------------
sys.Nprt                    = 5e3;                                          % Number of particles
sys.fun                     = @TwTaTm_aIRcoRoiRb;                           % model, Figure 5.1
sys.dt                      = dt;                                           % sampling time
sys.y                       = y;                                            % output
sys.u                       = u;                                            % input
sys.ahold                   = 'foh';                                        % first order hold
sys.sc                      = 1e8;                                          % scaling of thermal capacities
sys.H                       = [0 1 0];                                      % output matrix
Np                          = 9;                                            % Number of R,C parameters
Nx                          = 3;                                            % Number of states
Ny                          = 1;                                            % Number of outputs
Nu                          = 5;                                            % Number of inputs
sys.dim                     = [Np Nx Ny Nu];                                % Pass dimensions
%                             [ lb     theta     ub ]
Ro                          = [1e-4 ; 2.5e-2 ; 2e-1];
Ri                          = [1e-4 ; 2.5e-2 ; 2e-1];
Rm                          = [1e-4 ; 5e-3   ; 2e-1];
Rb                          = [1e-4 ; 5e-3   ; 2e-1];
Cw                          = [1e-4 ; 1e-1   ;  1e1];
Ca                          = [1e-4 ; 5e-3   ;  1e1];
Cm                          = [1e-4 ; 2e-1   ;  1e1];
aI                          = [1e-4 ; 4e-1   ;  1e0];
Rso                         = [1e-4 ; 3e-3   ; 2e-1];
gather                      = [Ro Ri Rm Rb Cw Ca Cm aI Rso];                % gather in one variable
m0                          = gather(2,:)';                                 % get parameters
sys.lb                      = gather(1,:)';                                 % get lower bounds
sys.ub                      = gather(3,:)';                                 % get upper bounds
sys.alpha                   = 0.05;                                         % 5% of parameter evolution authorized
m0std                       = [5e-3 5e-3 1e-3 1e-3 2e-2 ...                 % initial parameter std
                                    1e-3 3e-2 7e-2 5e-4]';                  % ...
FIM                         = diag([1e-6 1e-6 1e-6 1e-6 ...                 % Initial inverse Fisher information matrix
                                    1e-6 1e-6 1e-6 1e-4 1e-6]);             % ...
x0                          = [27 y(1) 27]';                                % initial states
P0                          = diag([0.5 0.25 0.5]);                         % initial state covariance
sys.P                       = repmat(P0,[1 1 sys.Nprt]);                    % replicate P0 Nprt times
sys.FIM                     = repmat(FIM, [1 1 sys.Nprt]);                  % replicate FIM Nprt times
sys.sig                     = diag([1e-5 1e-4 1e-4]);                       % scaling of incremental variance
sys.R                       = 1e-5;                                         % measurement noise covariance
sys.x                       = bsxfun(@plus,x0,bsxfun(@times, ...            % States particles
                                    sqrt(diag(P0)),randn(Nx,sys.Nprt)));    % ...
sys.m                       = bsxfun(@plus,m0,bsxfun(@times, ...            % Parameter particles
                                    m0std,randn(Np,sys.Nprt)));             % ...
out                         = Proposed(sys);                                % Proposed algorithm
xS                          = out.xS;                                       % get desired output from the structure
m                           = out.m;                                        % ...
