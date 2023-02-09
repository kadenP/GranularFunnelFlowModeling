%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bin insulation performance and sensitivity testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
Ins = Insulation();

%% Set parameters for full-scale storage bin
f = 1;
Ins.Hp = 7*f;   % m
Ins.bp = 2.25*f/Ins.Hp; 

% set simulation parameters
Ins.Tinf = 21;
Ins.T0 = 800;
Ins.thetaA = 0.5;
Ins.H = 0.1;

% set mesh parameters 
Ins.nrbar = 500;
Ins.nzbar = 100;

%% compute heat transfer properties
Ins.hInf = 10;   % (W/m2-K) heat transfer coefficient to ambient
Ins.hcw = 50;     % (W/m2-K) contact resistance coefficient
Ins.hcwA = 10;     % (W/m2-K) contact resistance coefficient
Ins.kp = 0.4;      % (W/mK) see Baumann and Zunft

% set insulation specifications for bin wall

Ins.nzbarW = 1000;
Ins.wallInsulation{1, 1} = 'tufcrete 47';
Ins.nrbarW{1} = 500;
Ins.wallInsulation{1, 2} = [2.25*f, 2.35*f];    % m
Ins.wallInsulation{1, 3} = 1.53;            % W/mK
Ins.wallInsulation{1, 4} = 2210;            % kg/m3
Ins.wallInsulation{1, 5} = 1175;            % J/kgK
Ins.wallInsulation{1, 6} = 100;             % number of lumped elements
Ins.wallInsulation{1, 7} = 25;              % W/m2K contact conductance
Ins.wallInsulation{1, 8} = 0;               % nondimensional IC

Ins.wallInsulation{2, 1} = 'skamolex';
Ins.nrbarW{2} = 500;
Ins.wallInsulation{2, 2} = [2.35*f, 2.55*f];        % m
Ins.wallInsulation{2, 3} = 0.09;                % W/mK
Ins.wallInsulation{2, 4} = 245;                 % kg/m3
Ins.wallInsulation{2, 5} = 840;                 % J/kgK
Ins.wallInsulation{2, 6} = 200;                 % number of lumped elements
Ins.wallInsulation{2, 7} = 100;                 % W/m2K contact conductance
Ins.wallInsulation{2, 8} = 0;                   % nondimensional IC

Ins.wallInsulation{3, 1} = 'elmtherm';
Ins.nrbarW{3} = 200;
Ins.wallInsulation{3, 2} = [2.55*f, 2.65*f];        % m
Ins.wallInsulation{3, 3} = 0.025;               % W/mK
Ins.wallInsulation{3, 4} = 270;                 % kg/m3
Ins.wallInsulation{3, 5} = 1005;                % J/kgK
Ins.wallInsulation{3, 6} = 100;                 % number of lumped elements
Ins.wallInsulation{3, 7} = 100;                 % W/m2K contact conductance
Ins.wallInsulation{3, 8} = 0;                   % nondimensional IC

Ins.wallInsulation{4, 1} = 'ss304';
Ins.nrbarW{4} = 50;
Ins.wallInsulation{4, 2} = [2.65*f, 2.65635*f];     % m
Ins.wallInsulation{4, 3} = 30;                  % W/mK
Ins.wallInsulation{4, 4} = 7700;                % kg/m3
Ins.wallInsulation{4, 5} = 500;                 % J/kgK
Ins.wallInsulation{4, 6} = 10;                  % number of lumped elements
Ins.wallInsulation{4, 7} = 1000;                % W/m2K contact conductance
Ins.wallInsulation{4, 8} = 0;                   % nondimensional IC

% % nutec block option

% Ins.wallInsulation{1, 1} = 'nutec';
% Ins.nrbarW{1} = 400;
% Ins.wallInsulation{1, 2} = [2.25*f, 2.6564*f];    % m
% Ins.wallInsulation{1, 3} = 0.22;              % W/mK
% Ins.wallInsulation{1, 4} = 160 ;              % kg/m3
% Ins.wallInsulation{1, 5} = 1130;              % J/kgK
% Ins.wallInsulation{1, 6} = 300;               % number of lumped elements
% Ins.wallInsulation{1, 7} = 25;                % W/m2K contact conductance
% Ins.wallInsulation{1, 8} = 0;                   % nondimensional IC
% 
% Ins.wallInsulation{2, 1} = 'ss304';
% Ins.nrbarW{2} = 50;
% Ins.wallInsulation{2, 2} = [2.6564*f, 2.6627*f];  % m
% Ins.wallInsulation{2, 3} = 30;                % W/mK
% Ins.wallInsulation{2, 4} = 7700;              % kg/m3
% Ins.wallInsulation{2, 5} = 500;               % J/kgK
% Ins.wallInsulation{2, 6} = 10;                % number of lumped elements
% Ins.wallInsulation{2, 7} = 1000;              % W/m2K contact conductance
% Ins.wallInsulation{2, 8} = 0;                   % nondimensional IC

% set insulation specifications for bin base
Ins.baseInsulation{1, 1} = 'particles';
Ins.baseInsulation{1, 2} = [0*f, 0.1*f];
Ins.baseInsulation{1, 3} = 0.4;             % W/mK
Ins.baseInsulation{1, 4} = 2000;            % kg/m3
Ins.baseInsulation{1, 5} = 1025.965;        % J/kgK
Ins.baseInsulation{1, 6} = 200;             % number of lumped elements
Ins.baseInsulation{1, 7} = 1e6;             % W/m2K contact conductance
Ins.baseInsulation{1, 8} = 0;               % nondimensional IC

Ins.baseInsulation{2, 1} = 'fondag';
Ins.baseInsulation{2, 2} = [0.1*f, (0.1+0.1905)*f];
Ins.baseInsulation{2, 3} = 1.75;          % W/mK
Ins.baseInsulation{2, 4} = 2210;          % kg/m3
Ins.baseInsulation{2, 5} = 1046.7;        % J/kgK
Ins.baseInsulation{2, 6} = 200;           % number of lumped elements
Ins.baseInsulation{2, 7} = 25;            % W/m2K contact conductance
Ins.baseInsulation{2, 8} = 0;             % nondimensional IC


% set insulation specifications for bin top
% new insulation configuration
Ins.roofInsulation{1, 1} = 'nutec';         % name
Ins.roofInsulation{1, 2} = [0, 0.3048*f];     % (m) dimensions/thickness [y1, y2]
Ins.roofInsulation{1, 3} = 0.4;             % (W/mK) thermal conductivity
Ins.roofInsulation{1, 4} = 64;              % (kg/m3) density
Ins.roofInsulation{1, 5} = 1130.4;          % (J/kgK) specific heat capacity

Ins.roofInsulation{2, 1} = 'elmtherm';
Ins.roofInsulation{2, 2} = [0.3048*f, (0.3048 + 0.0254)*f];
Ins.roofInsulation{2, 3} = 0.025;
Ins.roofInsulation{2, 4} = 270;
Ins.roofInsulation{2, 5} = 1005;

Ins.roofInsulation{3, 1} = 'steel ceiling';
Ins.roofInsulation{3, 2} = [(0.3048 + 0.0254)*f, (0.3048 + 0.0254 + 0.00635)*f];
Ins.roofInsulation{3, 3} = 30;
Ins.roofInsulation{3, 4} = 7700;
Ins.roofInsulation{3, 5} = 500;

%% Simulate simple step responses until s.s. is reached

% set timing parameters
Ins.df = Ins.t2Fo(60, 1);
Ins.FoEnd = Ins.t2Fo(3600*96, 1);
Ins.ztop = 0.8;
reInitObj(Ins);

% simulate wall step response
% [Twall, qWall] = simulateLWUS(Ins);
% [Tpw, qPW] = simulateLPW(Ins);
% [Tpw, qPW] = simulateLPSSW(Ins);
% [Twall, qWall, Tsswall, qssWall, e, r, t] = simulateLPSSWerror(Ins, 1);

% cyclic step response
tHigh = 16*3600;
tLow = 8*3600;
T0_ = 200;
dT_ = 200;
[Twall, qWall] = simulateCLWUS(Ins, tHigh, tLow, dT_, T0_);

% % simulate base step response
% [Tbase, qBase] = simulateLBUS(Ins);
% 
% % simulate roof step response
% [Troof, TwallA, Ta, qRadRoof, qRadWall] = simulateLRUS(Ins);

%% Simulate continuouse 1D wall model with particle domain

% % set timing parameters
% Ins.df = Ins.t2Fo(3600, 1);
% Ins.FoEnd = Ins.t2Fo(3600*48, 1);
% Ins.ztop = 0.8;
% reInitObj(Ins);
% 
% % simulate thermal decay in time
% buildWallCWP1D(Ins);
% computeCWm(Ins, 1);
% computeCWP1D(Ins, Ins.Fo);








