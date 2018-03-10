% Frequentist Model selection, Appendix.4 
% likelihood ratio between model M4 and M43, step 11 Table 4.6 
clear all;
close all;
clc
addpath(userpath, 'models');
addpath(userpath, 'models_ident');
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
Srgh                        = imp.data(:,55);                               % global solar radiation measured on a horizontal surface [W/m²]
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
Qhvac1                      = Qliving + Qchild  + Qbath;                    % Electric power injected in the south zone [W]
Tzone1_all                  = [Tliving Tbath Tchild Tcorridor]*Vspaces;     % south zone temperature, weighted mean according to ground surfaces [°C]
Tzone2_all                  = [Tkitchen Tdoorway Tparents]*Vspaces2;        % north zone temperature, weighted mean according to ground surfaces [°C]
%% ------------------------------------------------------------------------
% Structure
% -------------------------------------------------------------------------
sys.inputs                  = [Text(cut) Tzone2_all(cut) ...                % DIfferent inputs depending on the model
                               Srgh(cut) Qhvac1(cut)]';                     % ...
sys.inputs2                 = [Text(cut) Tzone2_all(cut) Sr_west(cut) ...   % ...
                               Sr_east(cut) Sr_south(cut) Qhvac1(cut)]';    % ...
sys.y                       = Tzone1_all(cut);                              % output
sys.ahold                   = 'foh';                                        % first order hold
sys.sc                      = 1e8;                                          % scaling of thermal capacities
sys.dt                      = dt;                                           % sampling time
sys.options                 = optimoptions(@fminunc,...                     % optimization options
                                'Display','iter-detailed',...               % ...
                                'Algorithm','quasi-newton',...              % ...
                                'SpecifyObjectiveGradient',true,...         % ...
                                'HessUpdate','bfgs',...                     % ...
                                'MaxIterations',1e3,...                     % ...
                                'MaxFunctionEvaluations',1e3,...            % ...
                                'OptimalityTolerance',1e-6,...              % ...
                                'StepTolerance',1e-6);                      % ...
Msub                        = feval(@M4_ident,sys);                         % smaller model
Np_sub                      = sum(Msub.fix);
Ns                          = length(sys.y);
[AIC_sub,BIC_sub]           = AIC_BIC(Msub.LL,Np_sub,Ns);

M                           = feval(@M43_ident,sys);                        % larger model, very long to converge
Np                          = sum(M.fix);
[AIC,BIC]                   = AIC_BIC(M.LL,Np,Ns);

pvalue                      = LLRatioTest(M.LL,Msub.LL,Np-Np_sub);          % pvalue LRT

% It is possible to evaluate a number of models by using a loop s.t.
% mdl = {@M21_ident @M22_ident @M23_ident ...}
% for i = 1:length(All)
%   M(i) = feval(mdl{i},sys)
% end
% It's a pseudo code and not a real implementation !!! 
