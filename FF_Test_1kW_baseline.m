%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 Kwth Small Bin Model
% Kaden Plewe
% 10/14/2022
% code used to generate predicted results for the small bin NSTTF tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
FF_test = FF();

%% Set parameters for full-scale storage bin
% prototype parameters for 8 hr charge, 10 hr hold, 6 hr discharge
FF_test.Hp = 0.85;   % m
% FF_test.bp = 2.25/FF_test.Hp;
FF_test.bp = 0.27/FF_test.Hp; 
mtot = 0.95*FF_test.Hp*pi*(FF_test.bp*FF_test.Hp)^2*FF_test.rhopPack;
Qcharge = 9*0.85/7; %mtot/(2*3600);
Qdischarge = 5*0.85/7; %mtot/(2*3600);
FF_test.Qp = Qdischarge/FF_test.rhopPack;
FF_test.QChp = Qcharge/FF_test.rhopPack;

% set simulation parameters
FF_test.Tinf = 21;
FF_test.T0 = 800;
FF_test.thetaA = 0.5;
FF_test.H = 0.1;
FF_test.h = 0.0011;
FF_test.a0 = 0.044;
FF_test.aOutlet = 0.044;
FF_test.a = 0.044;
FF_test.amc = 0.044;
FF_test.b = FF_test.bp;


FF_test.nrbar = 1000;
FF_test.nzbar0 = 1500;
FF_test.nzH0 = 2000;
FF_test.nzc0 = 250;
FF_test.nzhat = 15;
FF_test.nrhat = 25;
FF_test.nrtop = 100;

% FF_test.Uinf = 0.02; %0.0138;
FF_test.Uinfp = FF_test.Qp/(pi*(FF_test.a0*FF_test.Hp)^2);

FF_test.CapBase = 0.5e2;
FF_test.CapWall = 1e2;

FF_test.dz = 0.01;
FF_test.dr = 0.01;


% FF_test.mu = 1;
FF_test.dzmc = 0.002;
FF_test.drtop = 0.005;
FF_test.clim = 1e-7;     % only fourier coefficients > clim are used
FF_test.p = 4000;        % total number of beta values computed
FF_test.pf = 2500;       % range for beta values to be computed
FF_test.q = 4000;        % total number of eta values computed
FF_test.qf = 2500;       % range for eta values to be computed
FF_test.ni = 30;         % number of beta values used in computation of Cnm
FF_test.mi = 40;         % number of eta values used in computation of Cnm
climFD = 5e-7;           % only fourier coefficients > clim are used
pFD = 4000;              % total number of beta values computed
pfFD = 2500;             % range for beta values to be computed
qFD = 4000;              % total number of eta values computed
qfFD = 2500;             % range for eta values to be computed
niFD = 30;               % number of beta values used in computation of Cnm
miFD = 30;               % number of eta values used in computation of Cnm

FF_test.climH = 1e-7;    % only fourier coefficients > clim are used
FF_test.pH = 4000;       % total number of beta values computed
FF_test.pfH = 2500;      % range for beta values to be computed
FF_test.qH = 2500;       % total number of eta values computed
FF_test.qfH = 1000;      % range for eta values to be computed
FF_test.niH = 25;        % number of beta values used in computation of Cnm
FF_test.miH = 20;        % number of eta values used in computation of Cnm
FF_test.climFH = 5e-7;   % only fourier coefficients > clim are used
FF_test.pFH = 1000;      % total number of beta values computed
FF_test.pfFH = 100;      % range for beta values to be computed
FF_test.niFH = 25;       % number of beta values used in computation of Cnm

FF_test.climW = 1e-7;   % only fourier coefficients > clim are used
FF_test.pW = 500;      % total number of beta values computed
FF_test.pfW = 2000;     % range for beta values to be computed
FF_test.qW = 1000;      % total number of eta values computed
FF_test.qfW = 9000;     % range for eta values to be computed
FF_test.niW = 20;       % number of beta values used in computation of Cnm
FF_test.miW = 400;      % number of eta values used in computation of Cnm

FF_test.reInitObj;
%% compute heat transfer properties
FF_test.hInf = 20;   % (W/m2-K) heat transfer coefficient to ambient
FF_test.hcw = 10;     % (W/m2-K) contact resistance coefficient
FF_test.hcwA = 10;     % (W/m2-K) contact resistance coefficient
FF_test.kp = 0.4;      % (W/mK) see Baumann and Zunft

% set insulation specifications for bin wall
% new insulation configuration
FF_test.nzbarW = 1000;
FF_test.wallInsulation{1, 1} = 'tufcrete 47';
FF_test.nrbarW{1} = 100;
FF_test.wallInsulation{1, 2} = [0.27, 0.3145]; % m
FF_test.wallInsulation{1, 3} = 1.53;         % W/mK
FF_test.wallInsulation{1, 4} = 2210;         % kg/m3
FF_test.wallInsulation{1, 5} = 1175;         % J/kgK
FF_test.wallInsulation{2, 1} = 'skamolex';
FF_test.nrbarW{2} = 100;
FF_test.wallInsulation{2, 2} = [0.3945, 0.3399]; % m
FF_test.wallInsulation{2, 3} = 0.09;         % W/mK
FF_test.wallInsulation{2, 4} = 245;          % kg/m3
FF_test.wallInsulation{2, 5} = 840;          % J/kgK
FF_test.wallInsulation{3, 1} = 'elmtherm';
FF_test.nrbarW{3} = 50;
FF_test.wallInsulation{3, 2} = [0.3399, 0.3446];  % m
FF_test.wallInsulation{3, 3} = 0.025;        % W/mK
FF_test.wallInsulation{3, 4} = 270;          % kg/m3
FF_test.wallInsulation{3, 5} = 1005;         % J/kgK
FF_test.wallInsulation{4, 1} = 'ss304';
FF_test.nrbarW{4} = 10;
FF_test.wallInsulation{4, 2} = [0.3446, 0.35]; % m
FF_test.wallInsulation{4, 3} = 30;             % W/mK
FF_test.wallInsulation{4, 4} = 7700;           % kg/m3
FF_test.wallInsulation{4, 5} = 500;            % J/kgK


% set insulation specifications for bin base
% new insulation configuration
FF_test.baseInsulation{1, 1} = 'particles';
FF_test.baseInsulation{1, 2} = [0, 0.069];
FF_test.baseInsulation{1, 3} = 0.4;          % W/mK
FF_test.baseInsulation{1, 4} = 2000;          % kg/m3
FF_test.baseInsulation{1, 5} = 1025.965;        % J/kgK
FF_test.baseInsulation{2, 1} = 'fondag';
FF_test.baseInsulation{2, 2} = [0.069, (0.069+0.127)];
FF_test.baseInsulation{2, 3} = 1.75;          % W/mK
FF_test.baseInsulation{2, 4} = 2210;          % kg/m3
FF_test.baseInsulation{2, 5} = 1046.7;        % J/kgK
% FF_test.baseInsulation{3, 1} = 'ss304';
% FF_test.baseInsulation{3, 2} = [0.05+0.1905, 0.05+0.1905+0.0063];
% FF_test.baseInsulation{3, 3} = 30;          % W/mK
% FF_test.baseInsulation{3, 4} = 7700;        % kg/m3
% FF_test.baseInsulation{3, 5} = 500;         % J/kgK

% set insulation specifications for bin top
% new insulation configuration
FF_test.roofInsulation{1, 1} = 'nutec';         % name
FF_test.roofInsulation{1, 2} = [0, 0.3048]*FF_test.Hp/7;     % (m) dimensions/thickness [y1, y2]
FF_test.roofInsulation{1, 3} = 0.4;             % (W/mK) thermal conductivity
FF_test.roofInsulation{1, 4} = 64;              % (kg/m3) density
FF_test.roofInsulation{1, 5} = 1130.4;          % (J/kgK) specific heat capacity

FF_test.roofInsulation{2, 1} = 'elmtherm';
FF_test.roofInsulation{2, 2} = [0.3048, 0.3048 + 0.0254]*FF_test.Hp/7;
FF_test.roofInsulation{2, 3} = 0.025;
FF_test.roofInsulation{2, 4} = 270;
FF_test.roofInsulation{2, 5} = 1005;

FF_test.roofInsulation{3, 1} = 'steel ceiling';
FF_test.roofInsulation{3, 2} = [0.3048 + 0.0254, 0.3048 + 0.0254 + 0.00635]*FF_test.Hp/7;
FF_test.roofInsulation{3, 3} = 30;
FF_test.roofInsulation{3, 4} = 7700;
FF_test.roofInsulation{3, 5} = 500;

% compute prototype overall heat transfer coefficients to surroundings
computeUbase(FF_test); computeUwall(FF_test);
% initializeBaseSys(FF_test); initializeWallSys(FF_test);
FF_test.hp1 = 5;
FF_test.hp3 = 4;
FF_test.hp5D = 0.2;
FF_test.hp5H = 0.2;
FF_test.hp5C = 0.1;

reInitObj(FF_test);

%% match similarity parameters

% Reynolds number for momentum scaling (finding optimal viscosity)
FF_test.nu = FF_test.nup;
FF_test.Uinf = FF_test.nu*FF_test.Rep/FF_test.H;

% Prandlt number for energy scaling (finding optimal thermal diffusivity)
FF_test.alphaPacked = FF_test.nu/FF_test.Prp;

% thermal conductivity is relatively unaffected by particle selection at
% high temperatures, so keep equivalent for model and prototype
FF_test.k = FF_test.kp;
FF_test.cp = FF_test.cpp;

% reset model flow rates
FF_test.rhoPack = FF_test.k/(FF_test.alphaPacked*FF_test.cp);
FF_test.Q = FF_test.Uinf*pi*(FF_test.a0*FF_test.H)^2;
FF_test.QCh = Qcharge/Qdischarge*FF_test.Q;

% Biot numbers for boundary condition scaling
FF_test.h1 = FF_test.Bip1*FF_test.k/FF_test.H;
FF_test.h2 = FF_test.Bip2*FF_test.k/FF_test.H;
FF_test.h3 = FF_test.Bip3*FF_test.k/FF_test.H;
FF_test.h4 = FF_test.Bip4*FF_test.k/FF_test.H;
FF_test.h5D = FF_test.Bip5D*FF_test.k/FF_test.H;
FF_test.h5H = FF_test.Bip5H*FF_test.k/FF_test.H;
FF_test.h5C = FF_test.Bip5C*FF_test.k/FF_test.H;

reInitObj(FF_test);


%% load initial condition from previous cycle
% filepath = 'Cyclic Runs\Parametric Analysis\limiting cases\Test 6\Hold\theta_1.mat';
% thetaPrev = load(filepath, 'thetaSave');
% hold_start = thetaPrev.thetaSave;
% FF_test.thetaW = hold_start{end, 7};
% FF_test.rhoW = FF_test.thetaW;
% FF_test.thetaBase = hold_start{end, 8};
% hold_start = hold_start(end, :);


%% Simulate
testZ = FF_test.z;
testR = FF_test.r;

% set timing parameters
FF_test.df = FF_test.t2Fo(30, 1);
% FF_test.df = FF_test.t2Fo(60, 1);
FF_test.FoEnd = FF_test.t2Fo(60*(8+13+10), 1);

FF_test.thetaFolder = 'thetaTest';
FF_test.ls = 10;
FF_test.ztop = 0.1;
FF_test.reInitObj;
FF_test.deltaM = 0.005;

FF_test.thetaI = 0.9;

% dt = 60s

% for i = 1:8
%     FF_test.FoMode{1, i} = 'C';
% end
% for i = 8:21
%     FF_test.FoMode{1, i} = 'H';
% end
% for i = 21:length(FF_test.Fo)
%     FF_test.FoMode{1, i} = 'D';
% end

% dt = 30s

for i = 1:16
    FF_test.FoMode{1, i} = 'C';
end
for i = 16:42
    FF_test.FoMode{1, i} = 'H';
end
for i = 42:length(FF_test.Fo)
    FF_test.FoMode{1, i} = 'D';
end

simulateStorageTheta(FF_test, 0, 0);

%% simulate step response for wall and base models
t = 0:3600*24;
u = ones(size(t));
initializeWallSys(FF_test); initializeBaseSys(FF_test);
yWall = FF_test.computeWallSys(u, FF_test.t2Fo(t, 1));
yBase = FF_test.computeBaseSys(u, FF_test.t2Fo(t, 1));

figure('Units', 'normalized', 'color', 'white', ...
                                     'Position', [0 0 0.5 0.3]); hold on;
hold on;
plot(t/3600, FF_test.theta2T(yBase(:, 3)), '-k');
% title('Design 1 (Original)', 'interpreter', 'latex', 'FontSize', 14);
xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('$T$ ($^\circ$C)', 'interpreter', 'latex', 'FontSize', 14);






