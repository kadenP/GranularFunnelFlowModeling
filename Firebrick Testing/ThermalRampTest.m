%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integral solution class testing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
FF_test = FF();

%% Set parameters for full-scale storage bin
% prototype parameters for 8 hr charge, 10 hr hold, 6 hr discharge
FF_test.Hp = 7;   % m
% FF_test.bp = 2.25/FF_test.Hp;
FF_test.bp = 2.25/FF_test.Hp; 
mtot = 0.95*FF_test.Hp*pi*(FF_test.bp*FF_test.Hp)^2*FF_test.rhopPack;
Qcharge = 1.45; %mtot/(2*3600);
Qdischarge = 1.45; %mtot/(2*3600);
FF_test.Qp = Qdischarge/FF_test.rhopPack;
FF_test.QChp = Qcharge/FF_test.rhopPack;

% set simulation parameters
FF_test.Tinf = 21;
FF_test.T0 = 600;
FF_test.thetaA = 0.5;
FF_test.H = 0.1;
FF_test.h = 0.01; %0.0011;
FF_test.a0 = 0.044;
FF_test.aOutlet = 0.015;
FF_test.a = 0.044;
FF_test.amc = 0.044;
FF_test.b = FF_test.bp;
FF_test.dzhat = 0.0001; %0.001
FF_test.drbar = 0.00005;
FF_test.nrbar = 600;
% FF_test.drH = 0.0005;
FF_test.nrH = 600;
FF_test.drhat = 0.002;
% FF_test.Uinf = 0.02; %0.0138;
FF_test.Uinfp = FF_test.Qp/(pi*(FF_test.a0*FF_test.Hp)^2);

FF_test.CapBase = 0.5e2;
FF_test.CapWall = 1e2;

FF_test.dz = 0.01;
FF_test.dr = 0.01;

FF_test.dzbar = 0.0005;
FF_test.nzbar = 1000;
FF_test.dzH = 0.0005;
FF_test.nzH = 2000;
FF_test.dzc = 0.005;
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

FF_test.climH = 1e-6;    % only fourier coefficients > clim are used
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
FF_test.hInf = 10;   % (W/m2-K) heat transfer coefficient to ambient
FF_test.hcw = 50;     % (W/m2-K) contact resistance coefficient
FF_test.hcwA = 10;     % (W/m2-K) contact resistance coefficient
FF_test.kp = 0.4;      % (W/mK) see Baumann and Zunft

% set insulation specifications for bin wall
% new insulation configuration
FF_test.nzbarW = 1000;
FF_test.wallInsulation{1, 1} = 'tufcrete 47';
FF_test.nrbarW{1} = 100;
FF_test.wallInsulation{1, 2} = [2.25, 2.35]; % m
FF_test.wallInsulation{1, 3} = 1.53;         % W/mK
FF_test.wallInsulation{1, 4} = 2210;         % kg/m3
FF_test.wallInsulation{1, 5} = 1175;         % J/kgK
FF_test.wallInsulation{1, 6} = 100;          % number of lumped elements
FF_test.wallInsulation{1, 7} = 10;          % W/m2K contact conductance

FF_test.wallInsulation{2, 1} = 'skamolex';
FF_test.nrbarW{2} = 100;
FF_test.wallInsulation{2, 2} = [2.35, 2.55]; % m
FF_test.wallInsulation{2, 3} = 0.09;         % W/mK
FF_test.wallInsulation{2, 4} = 245;          % kg/m3
FF_test.wallInsulation{2, 5} = 840;          % J/kgK
FF_test.wallInsulation{2, 6} = 200;          % number of lumped elements
FF_test.wallInsulation{2, 7} = 1000;          % W/m2K contact conductance

FF_test.wallInsulation{3, 1} = 'elmtherm';
FF_test.nrbarW{3} = 50;
FF_test.wallInsulation{3, 2} = [2.55, 2.65];  % m
FF_test.wallInsulation{3, 3} = 0.025;        % W/mK
FF_test.wallInsulation{3, 4} = 270;          % kg/m3
FF_test.wallInsulation{3, 5} = 1005;         % J/kgK
FF_test.wallInsulation{3, 6} = 100;             % number of lumped elements
FF_test.wallInsulation{3, 7} = 50;           % W/m2K contact conductance

FF_test.wallInsulation{4, 1} = 'ss304';
FF_test.nrbarW{4} = 10;
FF_test.wallInsulation{4, 2} = [2.65, 2.65635]; % m
FF_test.wallInsulation{4, 3} = 30;             % W/mK
FF_test.wallInsulation{4, 4} = 7700;           % kg/m3
FF_test.wallInsulation{4, 5} = 500;            % J/kgK
FF_test.wallInsulation{4, 6} = 10;              % number of lumped elements
FF_test.wallInsulation{4, 7} = 1000;           % W/m2K contact conductance

% set insulation specifications for bin base
% new insulation configuration
FF_test.baseInsulation{1, 1} = 'particles';
FF_test.baseInsulation{1, 2} = [0, 0.1];
FF_test.baseInsulation{1, 3} = 0.4;          % W/mK
FF_test.baseInsulation{1, 4} = 2000;          % kg/m3
FF_test.baseInsulation{1, 5} = 1025.965;        % J/kgK
FF_test.baseInsulation{1, 6} = 50;          % number of lumped elements
FF_test.baseInsulation{1, 7} = 1;         % W/m2K contact conductance

FF_test.baseInsulation{2, 1} = 'fondag';
FF_test.baseInsulation{2, 2} = [0.1, 0.1+0.1905];
FF_test.baseInsulation{2, 3} = 1.75;          % W/mK
FF_test.baseInsulation{2, 4} = 2210;          % kg/m3
FF_test.baseInsulation{2, 5} = 1046.7;        % J/kgK
FF_test.baseInsulation{2, 6} = 100;           % number of lumped elements
FF_test.baseInsulation{2, 7} = 25;            % W/m2K contact conductance


% set insulation specifications for bin top
% new insulation configuration
FF_test.roofInsulation{1, 1} = 'nutec';         % name
FF_test.roofInsulation{1, 2} = [0, 0.3048];     % (m) dimensions/thickness [y1, y2]
FF_test.roofInsulation{1, 3} = 0.4;             % (W/mK) thermal conductivity
FF_test.roofInsulation{1, 4} = 64;              % (kg/m3) density
FF_test.roofInsulation{1, 5} = 1130.4;          % (J/kgK) specific heat capacity

FF_test.roofInsulation{2, 1} = 'steel ceiling';
FF_test.roofInsulation{2, 2} = [0.3048, 0.3048 + 0.00635];
FF_test.roofInsulation{2, 3} = 30;
FF_test.roofInsulation{2, 4} = 7700;
FF_test.roofInsulation{2, 5} = 500;

% compute prototype overall heat transfer coefficients to surroundings
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

%% Compute Charge and Discharge Times
tCharge = 0.98*mtot/Qcharge;
tHold = 4*3600;
tDischarge = 0.98*mtot/Qdischarge;


%% Simulate
testZ = FF_test.z;
testR = FF_test.r;

% set timing parameters
dt = 1200;
FF_test.df = FF_test.t2Fo(dt, 1);
FF_test.FoEnd = FF_test.t2Fo((tCharge + tHold + tDischarge), 1);

FF_test.thetaFolder = 'thetaTest';
FF_test.ls = 50;
FF_test.ztop = 0.05;
FF_test.reInitObj;
FF_test.deltaM = 0.005;

FF_test.tauW1 = FF_test.t2Fo(3600, 0);
FF_test.thetaI = 0.9;

for i = 1:ceil(tCharge/dt)
    FF_test.FoMode{1, i} = 'C';
end
for i = ceil(tCharge/dt):ceil((tCharge + tHold)/dt)
    FF_test.FoMode{1, i} = 'H';
end
for i = ceil((tCharge + tHold)/dt):length(FF_test.Fo)
    FF_test.FoMode{1, i} = 'D';
end

simulateStorageTheta(FF_test, 0, 0);
