classdef Insulation < handle
    % This class provides functions for an isolated technoeconomic analysis
    % of the bin insulation for the roof, base, and walls of a
    % flat-bottomed particle-based thermal energy storage bin. It uses
    % both lumped and spatially resolved models, enabling the user to test
    % the validity of the former.    
    properties
        % independent variables
        FoEnd = 1       % nondimensional simulation end time
        tEnd            % (s) simulation end time
        Fo              % nondimensional time (Fourier number) vector
        df = 0.1        % nondimensional time-step
        dt              % (s) time-step
        ztop = 0.9      % nondimensional height of particles
        FoNow = 0       % current iteration non-dimensional time
        zbar            % particle domain vertical mesh
        rbar            % particle domain radial mesh
        nzbar           % number of particle domain vertical mesh points
        nrbar           % number of particle domain radial mesh points
        zbarW           % z-dimension vectors for composite wall
        rbarW           % r-dimension vectors for composite wall
        nzbarW          % number of nodes used in zbarW
        nrbarW          % number of nodes used in rbarW
        dzbarW          % node spacing in zbarW
        drbarW          % node spacing in rbarW
        zbarB           % z-dimension vectors for composite base
        rbarB           % r-dimension vectors for composite base
        nzbarB          % number of nodes used in zbarB
        nrbarB          % number of nodes used in rbarB
        dzbarB          % node spacing in zbarB
        drbarB          % node spacing in rbarB
        zbarR           % z-dimension vectors for composite roof
        rbarR           % r-dimension vectors for composite roof
        nzbarR          % number of nodes used in zbarR
        nrbarR          % number of nodes used in rbarR
        dzbarR          % node spacing in zbarR
        drbarR          % node spacing in rbarR
        % prototype and model geometric and heat transfer parameters
        Hp = 7          % (m) height of prototype bin
        H = 1           % (m) height of stagnant region
        bp = 0.3214     % inner radius of prototype bin
        b = 0.3214      % outer nondimensional radius
        hInfp = 10      % (W/m2-K) prototype ambient heat transfer coefficient
        hInf            % (W/m2-K) model ambient heat transfer coefficient
        T0 = 800        % (°C) initial particle temperature
        Tinf = 20       % (°C) environmental ambient temperature
        % insulation layer information
        baseInsulation              % insulation info for base
        meshedBaseInsulation        % insulation info for discretized base
        zBL                         % (m) lumped vertical base mesh vector
        wallInsulation              % insulation info for wall
        meshedWallInsulation        % insulation info for discretized wall
        rWL                         % (m) lumped radial mesh vector
        roofInsulation              % insulation info for top of bin
        meshedRoofInsulation        % insulation info for discretized roof
        zRL                        % (m) lumped vertical roof mesh vector
        % insulation layer resistance and lumped capacitance variables
        Rwall           % (K/W) cell array storing resistance of wall layers
        Rcw             % (K/W) contact resistance between wall layers
        Rbase           % (K/W) cell array storing resistance of base layers
        Rcb             % (K/W) contact resistance between base layers
        Rroof           % (K/W) cell array storing resistance of top layers
        Rcr             % (K/W) contact resistance betweein roof layers
        ResP            % (K/W) convective resistance at particle top-air
        ResWA           % (K/W) cell array with resistance of exposed wall
        CapWall         % (J/K) cell array storing capacitance of wall layers
        CapBase         % (J/K) cell array storing capacitance of base layers
        CapRoof         % (J/K) cell array storing capacitance of top layers
        CapAir          % (J/K) thermal capacitance of internal air
        CapWA           % (J/K) cell array for capacitance of exposed wall
        % lumped model solution variables
        thetaA          % ambient temperature inside
        thetaWall       % state variables for wall
        thetaWA         % state variables for exposed wall
        thetaBase       % state variables for base
        thetaRoof       % state variables for top 
        % lumped model system matrices and system variables
        tauWall         % non-dimensional time constants for wall layers
        dTauwdr         % time constant sensitivity partial w.r.t thickness
        tauBase         % non-dimensional time constants for base layers
        tauTop          % non-dimensional time constants for top layers
        Awall           % state matrix for wall system
        dAwdr           % sensitivity state matrix w.r.t thickness
        Abase           % state matrix for base system
        Aroof           % state matrix for top system
        Bwall           % input matrix for wall system
        dBwdr           % input sensitivity matrix w.r.t thickness
        Bbase           % input matrix for base system
        Btop            % input matrix for top system
        Cwall           % state output matrix for wall system
        dCwdr           % state output sensitivity matrix for wall system
        Cbase           % state output matrix for base system
        Ctop            % state output matrix for top system
        Dwall           % input output matrix for wall system
        dDwdr           % input output sensitivity matrix for wall system
        Dbase           % input output matrix for base system
        Dtop            % input output matrix for top system
        wallSys         % wall state-space system
        baseSys         % base state-space system
        topSys          % top state-space system
        swr             % wall sensitivity index w.r.t thickness
        swrSys          % linear system for wall sensitivity index
        % enclosure radiation model variables and parameters
        epsilonP = 0.9      % emissivity of particle bed
        epsilonR = 0.8      % emissivity of exposed roof
        epsilonW = 0.8      % emissivity of exposed wall
        sigma = 5.67e-8     % (W/m2K4)
        cA = 1005           % (J/kgK) specific heat capacity of air
        tauR                % non-dimensional time constants in roof layers
        tauRA               % non-dimensional time constant roof-air
        tauAR               % non-dimensional time constant air-roof
        tauW                % non-dimensional time constants in wall layers
        tauWA               % non-dimensional time constant wall-air
        tauAW               % non-dimensional time constant air-wall
        tauAP               % non-dimensional time constant air-particles
        tauRinf             % non-dimensional time constant roof-ambient
        tauWinf             % non-dimensional time constant wall-ambient
        Ap                  % (m2) exposed area of particle bed
        Ar                  % (m2) exposed area of roof
        Aw                  % (m2) exposed area of wall
        Fpr                 % view factor between particles and roof
        Fpw                 % view factor between particles and wall
        Frw                 % view factor between roof and wall
        ha = 5      % (W/m2K) internal ambient heat transfer coefficient
        DAE         % set of all differential and algebraic equations
        AE          % set of algebraic equations
        DE          % set of differential equations
        Mass        % dae mass matrix
        qRadR       % (W) radiative heat flux into roof
        qRadW       % (W) radiative heat flux into exposed wall
        Tp          % (K) averaged temp on particle surface
        % spatially resolved model solution matrices
        thetaW          % composite wall solution/s 
        thetaW1D        % 1D composite wall solutions
        thetaWIC        % IC contribution for thetaW
        thetaWBC        % BC contribution for thetaW
        gW1             % current inner wall boundary condition
        gW1p            % averaged boundary condition with particles
        gW1a            % averaged boundary condition with air
        KWI             % boundary contribution kernel integrand
        KW              % boundary contribution kernel
        IBCW = []       % stored inner wall boundary conditions, f(z, t)
        rhoW            % composite wall initial condition  
        rhoW1D          % 1D composite wall initial condition
        % spatially resolved model coefficients and eigenvalues
        AWm             % r-BVP boundary condition matrix
        bWm             % r eigenfunction coefficients (IC contribution)
        bWBCm           % r eigenfunction coefficients (BC contribution)
        cWr             % lhs of reduced r-BVP boundary condition system
        etaW            % r-BVP eigenvalues reduced for IC contribution
        betaW           % z-BVP eigenvalues reduced for IC contribution
        etaWN           % r-BVP norms
        betaWN          % z-BVP norms
        etaWBC          % r-BVP eigenvalues reduced for BC contribution
        betaWBC         % z-BVP eigenvalues reduced for BC contribution
        betaWStatic     % static betaW values
        etaWStatic      % static etaW values
        AWmStatic       % boundary condition matrix for etaWStatic
        bWmStatic       % r-BVP coefficients for etaWStatic
        GWCnm           % wall coefficients for IC contribution
        CWm             % 1D wall coefficients
        GWBCCnm         % wall coefficients for BC contribution
        % other spatially resolved model parameters
        hcw = 1         % (W/m2-K) wall-particle boundary contact coefficient
        hcwA = 1        % (W/m2-K) wall-tank air convection coefficient
        Biw1            % biot number for inner-most composite layer
        Biw1A           % "" for tank air connection
        Biw2            % biot number for outer-most composite layer         
        % conduction fourier summation parameters
        climW = 5e-3     % only fourier coefficients > clim are used
        cGetW            % index array for obtaining eigenvalues
        qW = 8000        % total number of eta values computed
        qfW = 4000       % range for eta values to be computed
        miW = 2000       % number of eta values used in computation of Cnm
        % heat loss at domain boundaries
        qTopP           % (W/m2) heat loss from top of particle region
        qWall           % (W/m2) heat loss from wall composite layers
        qBase           % (W/m2) heat loss from base composite layers
        qLossW          % (W/m2) prescribed heat flux at bin wall
        qLossB          % (W/m2) prescribed heat flux at bin base
        qLossT          % (W/m2) prescribed heat flux at bin top
        qLossWp         % (kW) total wall heat loss computed in particles
        qLossBp         % (kW) total base heat loss computed in particles
        qLossTp         % (kW) total top heat loss computed in particles
        % particle parameters (for scaling)
        g = 9.80665                 % (m/s2)
        kp = 0.4                    % (W/mK) particle packed thermal 
                                    % conductivity
        rhopPack = 2000             % (kg/m3) particle packed bulk 
                                    % density
        cpp = 1025.965              % (J/kgK) average particle heat 
                                    % capacity 
        alphapPacked                % (m2/s) particle thermal 
                                    % diffusivity 
        k = 0.4                     % (W/mK) particle packed thermal 
                                    % conductivity
        rhoPack = 2000              % (kg/m3) particle packed bulk 
                                    % density
        cp = 1025.965               % (J/kgK) average particle heat 
                                    % capacity 
        alphaPacked                 % (m2/s) particle thermal 
                                    % diffusivity
        % conservation validation
        continuity      % cell-to-cell evaluation of 2D continuity equation
        energy          % cell-to-cell evaluation of 2D energy equation
        % data saving parameters
        ls = 10         % max storage size for time steps
        thetaFolder     % folder to save theta matrices in
    end
    
    methods
        function obj = Insulation()
            % scaling and property parameters
            obj.k = obj.kp;
            obj.cp = obj.cpp;
            obj.rhopPack = obj.rhoPack;
            obj.alphapPacked = obj.kp/(obj.rhopPack*obj.cpp); 
            obj.alphaPacked = obj.k/(obj.rhoPack*obj.cp);
            obj.hInf = obj.hInfp*(obj.Hp/obj.H);
            % radiation parameters
            obj.Ap = pi*(obj.bp*obj.Hp)^2;
            obj.Ar = obj.Ap;
            obj.Aw = 2*pi*obj.bp*obj.Hp*(obj.Hp*(1 - obj.ztop));
            r_ = obj.bp/(1 - obj.ztop);
            x = 2 + 1/r_^2;
            y = sqrt(x^2 - 4);
            obj.Fpr = 0.5*(x - y);
            obj.Fpw = 1 - obj.Fpr;
            obj.Frw = obj.Fpw;
            % timing parameters
            if isempty(obj.df) && ~isempty(obj.dt)
                obj.df = obj.t2Fo(obj.dt);
            elseif isempty(obj.dt) && ~isempty(obj.df)
                obj.dt = obj.Fo2t(obj.df);
            elseif ~isempty(obj.dt) && ~isempty(obj.df)
                if obj.t2Fo(obj.dt) ~= obj.df
                    warning('dt and df do not match, df used as default')
                    obj.dt = obj.Fo2t(obj.df);
                end
            else
                warning('time step not defined, default of 1s used');
                obj.dt = 1;
                obj.df = obj.t2Fo(obj.dt);
            end
            if isempty(obj.FoEnd) && ~isempty(obj.tEnd)
                obj.FoEnd = obj.t2Fo(obj.tEnd);
            elseif isempty(obj.tEnd) && ~isempty(obj.FoEnd)
                obj.tEnd = obj.Fo2t(obj.FoEnd);
            elseif ~isempty(obj.tEnd) && ~isempty(obj.FoEnd)
                if obj.t2Fo(obj.tEnd) ~= obj.FoEnd
                    warning('dt and df do not match, df used as default')
                    obj.dt = obj.Fo2t(obj.FoEnd);
                end
            else
                warning('end time not defined, default of 3600s used');
                obj.tEnd = 3600;
                obj.FoEnd = obj.t2Fo(obj.tEnd);
            end
            obj.Fo = linspace(0, obj.FoEnd, ceil(obj.FoEnd/obj.df))';                           
        end  
        function reInitObj(obj)
            % scaling and property parameters
            obj.k = obj.kp;
            obj.cp = obj.cpp;
            obj.rhopPack = obj.rhoPack;
            obj.alphapPacked = obj.kp/(obj.rhopPack*obj.cpp); 
            obj.alphaPacked = obj.k/(obj.rhoPack*obj.cp);
            obj.hInf = obj.hInfp*(obj.Hp/obj.H);
            % radiation parameters
            obj.Ap = pi*(obj.bp*obj.Hp)^2;
            obj.Ar = obj.Ap;
            obj.Aw = 2*pi*obj.bp*obj.Hp*(obj.Hp*(1 - obj.ztop));
            r_ = obj.bp/(1 - obj.ztop);
            x = 2 + 1/r_^2;
            y = sqrt(x^2 - 4);
            obj.Fpr = 0.5*(x - y);
            obj.Fpw = 1 - obj.Fpr;
            obj.Frw = obj.Fpw;
            % timing parameters
            if isempty(obj.df) && ~isempty(obj.dt)
                obj.df = obj.t2Fo(obj.dt);
            elseif isempty(obj.dt) && ~isempty(obj.df)
                obj.dt = obj.Fo2t(obj.df);
            elseif ~isempty(obj.dt) && ~isempty(obj.df)
                if obj.t2Fo(obj.dt) ~= obj.df
                    warning('dt and df do not match, df used as default')
                    obj.dt = obj.Fo2t(obj.df);
                end
            else
                warning('time step not defined, default of 1s used');
                obj.dt = 1;
                obj.df = obj.t2Fo(obj.dt);
            end
            if isempty(obj.FoEnd) && ~isempty(obj.tEnd)
                obj.FoEnd = obj.t2Fo(obj.tEnd);
            elseif isempty(obj.tEnd) && ~isempty(obj.FoEnd)
                obj.tEnd = obj.Fo2t(obj.FoEnd);
            elseif ~isempty(obj.tEnd) && ~isempty(obj.FoEnd)
                if obj.t2Fo(obj.tEnd) ~= obj.FoEnd
                    warning('dt and df do not match, df used as default')
                    obj.dt = obj.Fo2t(obj.FoEnd);
                end
            else
                warning('end time not defined, default of 3600s used');
                obj.tEnd = 3600;
                obj.FoEnd = obj.t2Fo(obj.tEnd);
            end
            obj.Fo = linspace(0, obj.FoEnd, ceil(obj.FoEnd/obj.df))';           
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % general simulation functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Twall, qWall] = simulateLWUS(obj, u)
            % simulates step response with a lumped composite wall model
            if nargin < 2, u = ones(size(obj.Fo)); end           
            initializeWallSys(obj);
            N = size(obj.meshedWallInsulation, 1);
            yWall = computeWallSys(obj, u, obj.Fo);
            obj.qWall = yWall(:, 1:N); qWall = obj.qWall;
            obj.qLossW = obj.qWall(:, 1);
            obj.thetaWall = [u, yWall(:, N+1:end)]; 
            Twall = obj.theta2T(obj.thetaWall);            
        end
        function [Twall, qWall] = simulateLPW(obj)
            % simulates step response with a lumped composite wall model
            u = zeros(size(obj.Fo));          
            initializeParticleWallSys(obj);
            N = size(obj.meshedWallInsulation, 1);
            yWall = computeWallSys(obj, u, obj.Fo);
            obj.qWall = yWall(:, 1:N); qWall = obj.qWall;
            obj.qLossW = obj.qWall(:, 1);
            obj.thetaWall = yWall(:, N+1:end); 
            Twall = obj.theta2T(obj.thetaWall);            
        end
        function [Twall, qWall] = simulateLPSSW(obj)
            % simulates step response with a lumped composite wall model
            if nargin < 2, u = ones(size(obj.Fo)); end           
            initializeParticleSys(obj);
            N = obj.nrbar;
            yP = computeParticleSys(obj, u, obj.Fo);            
            obj.qWall = yP(:, 1:N);
            obj.qLossW = obj.qWall(:, end);             
            [thetaW_, qW_] = computeSteadyStateWallTQ(obj, yP(:, end), ...
                                      obj.qLossW*(2*pi*obj.bp*obj.Hp^2));
            qWall = [obj.qWall, qW_];
            obj.thetaWall = [yP(:, N+1:end), thetaW_]; 
            Twall = obj.theta2T(obj.thetaWall);            
        end
        function [Tbase, qBase] = simulateLBUS(obj, u)
            % simulates step response with a lumped composite base model
            if nargin < 2, u = ones(size(obj.Fo)); end
            N = size(obj.baseInsulation, 1);
            initializeBaseSys(obj);
            yBase = computeBaseSys(obj, u, obj.Fo);
            obj.qBase = yBase(:, 1:N); qBase = obj.qBase;
            obj.qLossB = obj.qBase(:, 1);
            obj.thetaBase = yBase(:, N+1:end); 
            Tbase = obj.theta2T(obj.thetaBase);           
        end
        function [Twall, qWall, Tsswall, qssWall, e, r, t] = simulateLPSSWerror(obj, ShowPlot)
            % simulates the lumped transient wall model and steady state
            % wall model and computes the error according to the transient
            % baseline
            if nargin < 2, ShowPlot = 0; end
            [Twall, qWall] = simulateLPW(obj);
            data1 = animateDiscreteRadialLine(obj, 0.5, '', 0);
            [Tsswall, qssWall] = simulateLPSSW(obj);
            data2 = animateDiscreteRadialLine(obj, 0.5, '', 0);
            t = [data1{:, 1}]; r = data1{1, 2};
            [eT, eq] = plotSteadyStateModelError(obj, Twall, Tsswall, ...
                                           qWall, qssWall, t, r, ShowPlot);
            e = [eT, eq];
            if ShowPlot
                compareRadialTemps(obj, data1, data2, 0.5, ...
                        'ssRadialLineComp.gif', 'Transient Wall', ...
                        'Steady-State Wall');
            end            
        end
        function [Troof, Twall, Ta, qRadRoof, qRadWall] = simulateLRUS(obj)
            % simulates unit step response (unit step at top particle 
            % surface) with a lumped composite roof and radiation 
            % with the top particle surface and wall model
            N = size(obj.roofInsulation, 1); M = size(obj.wallInsulation, 1);
            obj.Tp = obj.theta2TK(1);
            initializeER(obj);
            Troof = zeros(length(obj.Fo), N);
            Twall = zeros(length(obj.Fo), M);
            Ta = zeros(size(obj.Fo));
            qRadRoof = zeros(size(obj.Fo));
            qRadWall = zeros(size(obj.Fo));
            Troof(1, :) = obj.theta2T(obj.thetaRoof);
            Twall(1, :) = obj.theta2T(obj.thetaWA);
            qRadRoof(1) = obj.qRadR;
            qRadWall(1) = obj.qRadW;
            Ta(1) = obj.theta2T(obj.thetaA);
            for i = 2:length(obj.Fo)
                computeER(obj, [0, obj.df]);
                Troof(i, :) = obj.theta2T(obj.thetaRoof);
                Twall(i, :) = obj.theta2T(obj.thetaWA);
                qRadRoof(i) = obj.qRadR;
                qRadWall(i) = obj.qRadW;
                Ta(i) = obj.theta2T(obj.thetaA);                                                
            end                        
        end
        function [Twall, qWall] = simulateCLWUS(obj, tHigh, tLow, dT_, T0_)
            % simulates the thermal response in the wall to a periodically
            % instantaneous charging and discharging particle boundary
            % cycle
            if nargin < 4, dT_ = 0; end
            if nargin < 5, T0_ = obj.T0; end
            u = zeros(size(obj.Fo));
            cycle = 1;
            for i = 1:length(obj.Fo)
                if mod(obj.Fo2t(obj.Fo(i), 1), tHigh + tLow) <= tHigh
                    u(i) = obj.T2theta(T0_ + dT_*(cycle - 1));
                else
                    u(i) = 0.5*obj.T2theta(T0_ + dT_*(cycle - 1));
                end 
                if i ~= 1
                    if mod(obj.Fo2t(obj.Fo(i), 1), tHigh + tLow) < ...
                            mod(obj.Fo2t(obj.Fo(i-1), 1), tHigh + tLow)
                        cycle = cycle + 1;
                    end
                end
            end
            initializeWallSys(obj);
            N = size(obj.meshedWallInsulation, 1);
            yWall = computeWallSys(obj, u, obj.Fo);
            obj.qWall = yWall(:, 1:N); qWall = obj.qWall;
            obj.qLossW = obj.qWall(:, 1);
            obj.thetaWall = [u, yWall(:, N+1:end)]; 
            Twall = obj.theta2T(obj.thetaWall);            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % storage bin top, base and wall RC model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeBaseSys(obj)
            % computes the linear RC system for the storage tank base
            setLumpedBaseMesh(obj);
            N = size(obj.meshedBaseInsulation, 1); 
            obj.Rbase = {}; obj.Rcb = {}; obj.CapBase = {}; 
            % insulation layer capacitance and resistance
            for i = 1:N
                obj.Rcb{i, 1} = 'contact';
                obj.Rbase{i, 1} = obj.meshedBaseInsulation{i, 1};
                obj.CapBase{i, 1} = obj.meshedBaseInsulation{i, 1};
                t = abs((obj.meshedBaseInsulation{i, 2}(2) - ...
                                 obj.meshedBaseInsulation{i, 2}(1)));
                k_ = obj.meshedBaseInsulation{i, 3};
                rho_ = obj.meshedBaseInsulation{i, 4};
                c_ = obj.meshedBaseInsulation{i, 5};
                hc_ = obj.meshedBaseInsulation{i, 6};
                obj.Rbase{i, 2} = t/(pi*(obj.b*obj.Hp)^2*k_);
                obj.Rcb{i, 2} = 1/(pi*(obj.bp*obj.Hp)^2*hc_);
                obj.CapBase{i, 2} = rho_*c_*t*pi*(obj.b*obj.Hp)^2;
            end
            % convective resistance
            obj.Rcb{N + 1, 1} = 'no contact';
            obj.Rcb{N + 1, 2} = 0;
            obj.Rbase{N + 1, 1} = 'convection';
            obj.Rbase{N + 1, 2} = 1/(pi*(obj.bp*obj.Hp)^2*obj.hInf);
            % time constants
            obj.tauBase = NaN*ones(N, 2);
            obj.tauBase(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                        [obj.CapBase{:, 2}].* ...
                        ([obj.Rbase{1:end-1, 2}] + [obj.Rcb{1:end-1, 2}]));
            obj.tauBase(:, 2) = obj.Hp^2./(obj.alphapPacked* ...
                        [obj.CapBase{:, 2}].* ...
                        ([obj.Rbase{2:end, 2}] + [obj.Rcb{2:end, 2}]));
            % state-space system
            obj.Abase = spdiags([-(obj.tauBase(:, 1) + obj.tauBase(:, 2)), ...
                        obj.tauBase(:, 1), obj.tauBase(:, 2)], ...
                        [0, 1, -1], N, N)';
            obj.Bbase = zeros(N, 1); obj.Bbase(1) = obj.tauBase(1, 1);
            obj.Cbase = [zeros(N); eye(N)]; 
            for i = 1:N
                if i == N
                    obj.Cbase(i, i) = -(obj.T0 - obj.Tinf)/ ...
                      ((obj.Rbase{i, 2} + obj.Rcb{i, 2})*pi*(obj.bp*obj.Hp)^2);
                    obj.Cbase(i, i-1) = (obj.T0 - obj.Tinf)/ ...
                        ((obj.Rbase{i, 2} + obj.Rcb{i, 2})*pi*(obj.bp*obj.Hp)^2);
                else
                    obj.Cbase(i, i) = (obj.T0 - obj.Tinf)/ ...
                        ((obj.Rbase{i, 2} + obj.Rcb{i, 2})*pi*(obj.bp*obj.Hp)^2);
                    obj.Cbase(i, i+1) = -(obj.T0 - obj.Tinf)/ ...
                        ((obj.Rbase{i, 2} + obj.Rcb{i, 2})*pi*(obj.bp*obj.Hp)^2);
                end
            end
            obj.Dbase = zeros(2*N, 1);
            obj.baseSys = ...
                ss(full(obj.Abase), obj.Bbase, obj.Cbase, obj.Dbase); 
        end
        function initializeWallSys(obj, sensitivity)
            % computes the linear RC system for the storage tank base
            if nargin < 2, sensitivity = 0; end
            setLumpedWallMesh(obj);
            N = size(obj.meshedWallInsulation, 1); 
            obj.Rwall = {}; obj.Rcw = {}; obj.CapWall = {}; 
            w_ = zeros(N, 1);
            % insulation layer capacitance and resistance
            for i = 1:N
                obj.Rcw{i, 1} = 'contact';
                obj.Rwall{i, 1} = obj.meshedWallInsulation{i, 1};
                obj.CapWall{i, 1} = obj.meshedWallInsulation{i, 1};
                r1 = obj.meshedWallInsulation{i, 2}(1);
                r2 = obj.meshedWallInsulation{i, 2}(2);
                w_(i) = r2 - r1;
                k_ = obj.meshedWallInsulation{i, 3};
                rho_ = obj.meshedWallInsulation{i, 4};
                c_ = obj.meshedWallInsulation{i, 5};
                hc_ = obj.meshedWallInsulation{i, 6};
                obj.Rwall{i, 2} = log(r2/r1)/(2*pi*obj.Hp*k_);
                obj.Rcw{i, 2} = 1/(2*pi*r1*obj.Hp*hc_);
                obj.CapWall{i, 2} = rho_*c_*obj.Hp*pi*(r2^2 - r1^2);
            end
            % convective resistance
            obj.Rcw{N + 1, 1} = 'no contact';
            obj.Rcw{N + 1, 2} = 0;
            obj.Rwall{N + 1, 1} = 'ambient convection';
            obj.Rwall{N + 1, 2} = 1/(2*pi*r2*obj.Hp*obj.hInf);
            % time constants
            obj.tauWall = NaN*ones(N, 2);
            obj.tauWall(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                         [obj.CapWall{:, 2}].* ...
                         ([obj.Rwall{1:end-1, 2}] + [obj.Rcw{1:end-1, 2}]));
            obj.tauWall(:, 2) = obj.Hp^2./(obj.alphapPacked* ...
                         [obj.CapWall{:, 2}].* ...
                         ([obj.Rwall{2:end, 2}] + [obj.Rcw{2:end, 2}]));
            % state-space system
            obj.Awall = spdiags([-(obj.tauWall(:, 1) + obj.tauWall(:, 2)), ...
                        obj.tauWall(:, 1), obj.tauWall(:, 2)], ...
                        [0, 1, -1], N, N)';
            obj.Bwall = zeros(N, 1); obj.Bwall(1) = obj.tauWall(1, 1);
            obj.Cwall = [zeros(N); eye(N)];
            for i = 1:N
                if i == N
                    obj.Cwall(i, i) = -(obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                    obj.Cwall(i, i-1) = (obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                else
                    obj.Cwall(i, i) = (obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                    obj.Cwall(i, i+1) = (obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                end
            end
            obj.Dwall = zeros(2*N, 1);
            obj.wallSys = ...
                ss(full(obj.Awall), obj.Bwall, obj.Cwall, obj.Dwall);
            % construct sensitivity model
            if sensitivity
                % compute the set of non-dimensional time constant partials
                obj.dTauwdr = zeros(N, N+1);
                for i = 1:N
                    for j = 1:i+1
                        if j == N + 1
                            obj.dTauwdr(i, j) = 0;
                        else
                            ri = obj.meshedWallInsulation{i, 2}(1);
                            rip1 = obj.meshedWallInsulation{i, 2}(2);
                            rj = obj.meshedWallInsulation{j, 2}(1);
                            rjp1 = obj.meshedWallInsulation{j, 2}(2);
                            wSumi = obj.b*obj.Hp + sum(w_(1:i));
                            wSumj = obj.b*obj.Hp + sum(w_(1:j));
                            if i == 1
                                wSumim1 = obj.b*obj.Hp;
                            else
                                wSumim1 = obj.b*obj.Hp + sum(w_(1:i-1));
                            end
                            if j == 1
                                wSumjm1 = obj.b*obj.Hp;
                            else
                                wSumjm1 = obj.b*obj.Hp + sum(w_(1:j-1));
                            end
                            kj = obj.meshedWallInsulation{j, 3};
                            rhoi = obj.meshedWallInsulation{i, 4};
                            ci = obj.meshedWallInsulation{i, 5};
                            if i == j
                                obj.dTauwdr(i, j) = 2*obj.Hp^2*kj/ ...
                                    (rhoi*ci*obj.alphapPacked) ...
                                    *-2*wSumi*(wSumi^2 - wSumim1^2)^-2 ...
                                    *(log(wSumj/wSumjm1))^-1 - 1/wSumj ...
                                    *(log(wSumj/wSumjm1))^-2*(wSumi^2 - wSumim1^2)^-1;
                            else
                               obj.dTauwdr(i, j) = 2*obj.Hp^2*kj/ ...
                                    (rhoi*ci*obj.alphapPacked) ...
                                    *-2*wSumi*(wSumi^2 - wSumim1^2)^-2 ...
                                    *(log(wSumj/wSumjm1))^-1; 
                            end 
                        end
                    end
                end
                % compute sensitivity matrices
                obj.dAwdr = cell(N, 1);
                obj.dBwdr = cell(N, 1);
                obj.dCwdr = cell(2*N, N);
                obj.dDwdr = cell(2*N, 1);
                for i = 1:N
                    obj.dAwdr{i} = zeros(N);
                    obj.dBwdr{i} = zeros(N, 1);
                    obj.dCwdr{i} = zeros(2*N, N);
                    obj.dDwdr{i} = zeros(2*N, 1);
                    for j = i:N
                        kj = obj.meshedWallInsulation{j, 3};
                        rj = obj.meshedWallInsulation{j, 2}(1);
                        rjp1 = obj.meshedWallInsulation{j, 2}(2);
                        if j == 1
                            obj.dAwdr{i}(j, j) = -(obj.dTauwdr(j, j) + ...
                                               obj.dTauwdr(j, j+1));
                            obj.dAwdr{i}(j, j+1) = obj.dTauwdr(j, j+1);
                            obj.dBwdr{i}(j) = obj.dTauwdr(j, j);
                            obj.dCwdr{i}(i, i) = kj*(obj.T0 - obj.Tinf) ...
                               /(rjp1*log(rjp1/rj)^2*obj.b*obj.Hp);
                            obj.dDwdr{i}(i) = -kj*(obj.T0 - obj.Tinf) ...
                               /(rjp1*log(rjp1/rj)^2*obj.b*obj.Hp);
                        elseif j == N
                            obj.dAwdr{i}(j, j-1) = obj.dTauwdr(j, j);
                            obj.dAwdr{i}(j, j) = -(obj.dTauwdr(j, j) + ...
                                               obj.dTauwdr(j, j+1)); 
                            obj.dCwdr{i}(j, j) = kj*(obj.T0 - obj.Tinf) ...
                               /(rjp1*log(rjp1/rj)^2*obj.b*obj.Hp);
                            obj.dCwdr{i}(j, j-1) = -kj*(obj.T0 - obj.Tinf) ...
                               /(rjp1*log(rjp1/rj)^2*obj.b*obj.Hp);
                        else
                            obj.dAwdr{i}(j, j-1) = obj.dTauwdr(j, j);
                            obj.dAwdr{i}(j, j) = -(obj.dTauwdr(j, j) + ...
                                               obj.dTauwdr(j, j+1));                        
                            obj.dAwdr{i}(j, j+1) = obj.dTauwdr(j, j+1); 
                            obj.dCwdr{i}(j, j) = kj*(obj.T0 - obj.Tinf) ...
                               /(rjp1*log(rjp1/rj)^2*obj.b*obj.Hp);
                            obj.dCwdr{i}(j, j-1) = -kj*(obj.T0 - obj.Tinf) ...
                               /(rjp1*log(rjp1/rj)^2*obj.b*obj.Hp);
                        end
                    end
                    A_ = full(obj.Awall);
                    B_ = [obj.dBwdr{i}, obj.dAwdr{i}];
                    C_ = obj.Cwall;
                    D_ = [obj.dDwdr{i}, obj.dCwdr{i}];
                    obj.swrSys{i} = ss(A_, B_, C_, D_);
                    obj.swr{i} = zeros(N, 1);
                end
            end
        end
        function y = computeBaseSys(obj, u, Fo_)
            % computes the heat flux leaving the
            if nargin < 2, u = ones(2, 1); end
            if nargin < 3, Fo_ = linspace(0, obj.df, length(u)); end
            y = lsim(obj.baseSys, u, Fo_, obj.thetaBase);
            obj.qLossB = y(end, 1); obj.thetaBase = y(end, 2:end);
        end
        function y = computeWallSys(obj, u, Fo_)
            % computes the heat flux leaving the wall
            if nargin < 2, u = ones(2, 1); end
            if nargin < 3, Fo_ = linspace(0, obj.df, length(u)); end
            N = size(obj.meshedWallInsulation, 1);
            y = lsim(obj.wallSys, u, Fo_, obj.thetaWall);
            obj.qLossW = y(end, 1); obj.thetaWall = y(end, N+1:end);
        end
        function y = computeWallSensitivityR(obj, u, Fo_)
            % computes the time propogation of the sensitivity w.r.t the
            % thickness of wall layers
            initializeWallSys(obj, 1);
            N = size(obj.wallInsulation, 1);
            yW_ = computeWallSys(obj, u, Fo_);
            x = yW_(:, N+1:end);
            u_ = [u; x'];
            y = cell(2*N, 1);
            for i = 1:N
                y{i} = lsim(obj.swrSys{i}, u_, Fo_, obj.swr{i});
                y{i}(:, 1:N) = y{i}(:, 1:N)/obj.H;
            end
        end
        function y = computeSteadyWallSensitivityR(obj, xss, uss)
            % computes the steady state sensitivity metric for a given
            % steady state and steady state input
            initializeWallSys(obj, 1);
            N = size(obj.wallInsulation, 1);
            y = cell(N, 1);
            for i = 1:N
               spss = -obj.Awall\([obj.dAwdr{i}, obj.dBwdr{i}]*[xss; uss]);
               y{i} = obj.Cwall*spss + [obj.dCwdr{i}, obj.dDwdr{i}]*[xss; uss];                
            end
        end
        function [y, qss, qppss] = computeSteadyWallFluxSensitivityR(obj)
            % computes the steady state sensitivity of the wall heat loss
            % to each of the wall layer thicknesses
            N = size(obj.wallInsulation, 1);
            Rtot = sum(cell2mat(obj.Rwall(:, 2)));
            qss = (obj.T0 - obj.Tinf)/Rtot;
            y = zeros(N, 1); qppss = zeros(N, 1);
            w_ = zeros(N, 1); wSumi = zeros(N, 1); wSumim1 = zeros(N, 1);
            for i = 1:N
                r1 = obj.wallInsulation{i, 2}(1);
                r2 = obj.wallInsulation{i, 2}(2);
                w_(i) = r2 - r1;
                qppss(i) = qss/(2*pi*r1*obj.Hp);
                wSumi(i) = obj.b*obj.Hp + sum(w_(1:i));
                if i == 1
                    wSumim1(i) = obj.b*obj.Hp;
                else
                    wSumim1(i) = obj.b*obj.Hp + sum(w_(1:i-1));
                end                
            end
            for i = 1:N
                for j = i:N
                    ki = obj.wallInsulation{i, 3};
                    if i == j
                        dRtotdwi = (wSumi(i)*2*pi*obj.Hp*ki)^(-1);                        
                    else
                        dRtotdwi = dRtotdwi + ...
                            (1/wSumi(i) - 1/wSumim1(i))/(2*pi*obj.Hp*ki);                        
                    end
                end
               y(i) = -Rtot^-2*(obj.T0 - obj.Tinf)*dRtotdwi;                 
            end
        end
        function setLumpedWallMesh(obj)
            % adds discrete lumped elements for the wall model to increase
            % spatial resolution
            N = size(obj.wallInsulation, 1);
            obj.meshedWallInsulation = {};
            p = 1;
            for i = 1:N
                M = obj.wallInsulation{i, 6};
                rmin = obj.wallInsulation{i, 2}(1);
                rmax = obj.wallInsulation{i, 2}(2);
                r_ = linspace(rmin, rmax, M+1);
                for j = 1:M
                    obj.meshedWallInsulation{p, 1} = ...
                        obj.wallInsulation{i, 1};
                    obj.meshedWallInsulation{p, 2} = [r_(j), r_(j+1)];
                    obj.meshedWallInsulation{p, 3} = ...
                        obj.wallInsulation{i, 3};
                    obj.meshedWallInsulation{p, 4} = ...
                        obj.wallInsulation{i, 4};
                    obj.meshedWallInsulation{p, 5} = ...
                        obj.wallInsulation{i, 5};
                    if j == 1
                        obj.meshedWallInsulation{p, 6} = ...
                        obj.wallInsulation{i, 7};
                    else
                        obj.meshedWallInsulation{p, 6} = 1e16;
                    end
                    obj.thetaWall(p) = obj.wallInsulation{i, 8};
                    p = p+1;
                end             
            end 
            obj.rWL = unique(cell2mat(obj.meshedWallInsulation(:, 2)), 'sorted');
        end
        function setLumpedBaseMesh(obj)
            % adds discrete lumped elements for the wall model to increase
            % spatial resolution
            N = size(obj.baseInsulation, 1);
            obj.meshedBaseInsulation = {};
            p = 1;
            for i = 1:N
                M = obj.baseInsulation{i, 6};
                zmin = obj.baseInsulation{i, 2}(1);
                zmax = obj.baseInsulation{i, 2}(2);
                z_ = linspace(zmin, zmax, M+1);
                for j = 1:M
                    obj.meshedBaseInsulation{p, 1} = ...
                        obj.baseInsulation{i, 1};
                    obj.meshedBaseInsulation{p, 2} = [z_(j), z_(j+1)];
                    obj.meshedBaseInsulation{p, 3} = ...
                        obj.baseInsulation{i, 3};
                    obj.meshedBaseInsulation{p, 4} = ...
                        obj.baseInsulation{i, 4};
                    obj.meshedBaseInsulation{p, 5} = ...
                        obj.baseInsulation{i, 5};
                    if j == 1
                        obj.meshedBaseInsulation{p, 6} = ...
                            obj.baseInsulation{i, 7};
                    else
                        obj.meshedBaseInsulation{p, 6} = 1e16;
                    end
                    obj.thetaBase(p) = obj.baseInsulation{i, 8};
                    p = p+1;
                end             
            end 
            obj.zBL = unique(cell2mat(obj.meshedBaseInsulation(:, 2)), 'sorted');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insulation model with transient particle domain
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeParticleWallSys(obj)
            % initializes the RC system for the particles and wall
            % insulation with an adiabatic centerline condition            
            setLumpedParticleWallMesh(obj);
            N = size(obj.meshedWallInsulation, 1); 
            obj.Rwall = {}; obj.Rcw = {}; obj.CapWall = {}; 
            w_ = zeros(N, 1);
            % insulation layer capacitance and resistance
            obj.Rcw{1, 1} = 'contact';
            obj.Rwall{1, 1} = obj.meshedWallInsulation{1, 1};
            obj.CapWall{1, 1} = obj.meshedWallInsulation{1, 1};
            r1 = obj.meshedWallInsulation{1, 2}(1);
            r2 = obj.meshedWallInsulation{1, 2}(2);
            w_(1) = r2 - r1;
            rho_ = obj.meshedWallInsulation{1, 4};
            c_ = obj.meshedWallInsulation{1, 5};
            obj.Rwall{1, 2} = NaN;
            obj.Rcw{1, 2} = NaN;
            obj.CapWall{1, 2} = rho_*c_*obj.Hp*pi*(r2^2 - r1^2);
            for i = 2:N
                obj.Rcw{i, 1} = 'contact';
                obj.Rwall{i, 1} = obj.meshedWallInsulation{i, 1};
                obj.CapWall{i, 1} = obj.meshedWallInsulation{i, 1};
                r1 = obj.meshedWallInsulation{i, 2}(1);
                r2 = obj.meshedWallInsulation{i, 2}(2);
                w_(i) = r2 - r1;
                k_ = obj.meshedWallInsulation{i, 3};
                rho_ = obj.meshedWallInsulation{i, 4};
                c_ = obj.meshedWallInsulation{i, 5};
                hc_ = obj.meshedWallInsulation{i, 6};
                obj.Rwall{i, 2} = log(r2/r1)/(2*pi*obj.Hp*k_);
                obj.Rcw{i, 2} = 1/(2*pi*r1*obj.Hp*hc_);
                obj.CapWall{i, 2} = rho_*c_*obj.Hp*pi*(r2^2 - r1^2);
            end
            % convective resistance
            obj.Rcw{N + 1, 1} = 'no contact';
            obj.Rcw{N + 1, 2} = 0;
            obj.Rwall{N + 1, 1} = 'ambient convection';
            obj.Rwall{N + 1, 2} = 1/(2*pi*r2*obj.Hp*obj.hInf);
            % time constants
            obj.tauWall = NaN*ones(N-1, 2);
            obj.tauWall(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                         [obj.CapWall{2:end, 2}].* ...
                         ([obj.Rwall{2:end-1, 2}] + [obj.Rcw{2:end-1, 2}]));
            obj.tauWall(:, 2) = obj.Hp^2./(obj.alphapPacked* ...
                         [obj.CapWall{2:end, 2}].* ...
                         ([obj.Rwall{3:end, 2}] + [obj.Rcw{3:end, 2}]));
            % state-space system
            obj.Awall = spdiags([-(obj.tauWall(:, 1) + obj.tauWall(:, 2)), ...
                        obj.tauWall(:, 1), obj.tauWall(:, 2)], ...
                        [0, 1, -1], N-1, N-1)';
            obj.Awall = [zeros(N, 1), [zeros(1, N-1); obj.Awall]];
            obj.Awall(2, 1) = obj.tauWall(1, 1);
            obj.Awall(1, :) = obj.Awall(2, :);
            obj.Bwall = zeros(N, 1);
            obj.Cwall = [zeros(N); eye(N)];
            for i = 1:N
                if i == 1
                    % no entry
                elseif i == N
                    obj.Cwall(i, i) = -(obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                    obj.Cwall(i, i-1) = (obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                else
                    obj.Cwall(i, i) = (obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                    obj.Cwall(i, i+1) = -(obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                end
            end
            obj.Dwall = zeros(2*N, 1);
            obj.wallSys = ...
                ss(full(obj.Awall), obj.Bwall, obj.Cwall, obj.Dwall);                        
        end
        
        function initializeParticleSys(obj)
            % initializes dynamic system for the particle domain with a
            % steady state wall insulation model
            setLumpedParticleWallMesh(obj);
            N = size(obj.meshedWallInsulation, 1);
            Np = obj.nrbar;
            obj.Rwall = {}; obj.Rcw = {}; obj.CapWall = {}; 
            w_ = zeros(N, 1);
            % insulation layer capacitance and resistance
            obj.Rcw{1, 1} = 'contact';
            obj.Rwall{1, 1} = obj.meshedWallInsulation{1, 1};
            obj.CapWall{1, 1} = obj.meshedWallInsulation{1, 1};
            r1 = obj.meshedWallInsulation{1, 2}(1);
            r2 = obj.meshedWallInsulation{1, 2}(2);
            w_(1) = r2 - r1;
            rho_ = obj.meshedWallInsulation{1, 4};
            c_ = obj.meshedWallInsulation{1, 5};
            obj.Rwall{1, 2} = NaN;
            obj.Rcw{1, 2} = NaN;
            obj.CapWall{1, 2} = rho_*c_*obj.Hp*pi*(r2^2 - r1^2);
            for i = 2:N
                obj.Rcw{i, 1} = 'contact';
                obj.Rwall{i, 1} = obj.meshedWallInsulation{i, 1};
                obj.CapWall{i, 1} = obj.meshedWallInsulation{i, 1};
                r1 = obj.meshedWallInsulation{i, 2}(1);
                r2 = obj.meshedWallInsulation{i, 2}(2);
                w_(i) = r2 - r1;
                k_ = obj.meshedWallInsulation{i, 3};
                rho_ = obj.meshedWallInsulation{i, 4};
                c_ = obj.meshedWallInsulation{i, 5};
                hc_ = obj.meshedWallInsulation{i, 6};
                obj.Rwall{i, 2} = log(r2/r1)/(2*pi*obj.Hp*k_);
                obj.Rcw{i, 2} = 1/(2*pi*r1*obj.Hp*hc_);
                obj.CapWall{i, 2} = rho_*c_*obj.Hp*pi*(r2^2 - r1^2);
            end
            % convective resistance
            obj.Rcw{N + 1, 1} = 'no contact';
            obj.Rcw{N + 1, 2} = 0;
            obj.Rwall{N + 1, 1} = 'ambient convection';
            obj.Rwall{N + 1, 2} = 1/(2*pi*r2*obj.Hp*obj.hInf);
            % construct set of particle data and ss overall convection model
            Rp_ = obj.Rwall(1:Np, :);
            Cp_ = obj.CapWall(1:Np, :);
            Rtot_ = sum([obj.Rwall{Np+1:end, 2}] + [obj.Rcw{Np+1:end, 2}]);
            Rp_{Np+1, 1} = 'overall convection coefficient';
            Rp_{Np+1, 2} = Rtot_;                       
            % time constants
            obj.tauWall = NaN*ones(Np-1, 2);
            obj.tauWall(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                         [Cp_{2:end, 2}].*[Rp_{2:end-1, 2}]);
            obj.tauWall(:, 2) = obj.Hp^2./(obj.alphapPacked* ...
                         [Cp_{2:end, 2}].*[Rp_{3:end, 2}]);
            % state-space system
            obj.Awall = spdiags([-(obj.tauWall(:, 1) + obj.tauWall(:, 2)), ...
                        obj.tauWall(:, 1), obj.tauWall(:, 2)], ...
                        [0, 1, -1], Np-1, Np-1)';
            obj.Awall = [zeros(Np, 1), [zeros(1, Np-1); obj.Awall]];
            obj.Awall(2, 1) = obj.tauWall(1, 1);
            obj.Awall(1, :) = obj.Awall(2, :);
            obj.Bwall = zeros(Np, 1);
            obj.Cwall = [zeros(Np); eye(Np)];
            for i = 1:Np
                if i == 1
                    % no entry
                elseif i == Np
                    obj.Cwall(i, i) = -(obj.T0 - obj.Tinf)/(Rp_{i, 2} ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                    obj.Cwall(i, i-1) = (obj.T0 - obj.Tinf)/(Rp_{i, 2} ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                else
                    obj.Cwall(i, i) = (obj.T0 - obj.Tinf)/(Rp_{i, 2} ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                    obj.Cwall(i, i+1) = -(obj.T0 - obj.Tinf)/(Rp_{i, 2} ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                end
            end
            obj.Dwall = zeros(2*Np, 1);
            obj.wallSys = ...
                ss(full(obj.Awall), obj.Bwall, obj.Cwall, obj.Dwall);           
        end
        function y = computeParticleSys(obj, u, Fo_)
            % computes the heat flux leaving the wall
            if nargin < 2, u = ones(2, 1); end
            if nargin < 3, Fo_ = linspace(0, obj.df, length(u)); end
            N = obj.nrbar;
            y = lsim(obj.wallSys, u, Fo_, obj.thetaWall(1:N));
            obj.qLossW = y(end, 1); obj.thetaWall = y(end, N+1:end);
        end
        function [thetaW_, qWall_] = computeSteadyStateWallTQ(obj, thetaP, qW)
            % computes the algebraic steady-state wall temps from the
            % transient particle temperatures
            N = size(obj.meshedWallInsulation, 1) - obj.nrbar;
            Np = obj.nrbar;
            thetaW_ = zeros(length(qW), N);
            qWall_ = zeros(length(qW), N);
            thetaW_(2:end, 1) = thetaP(2:end) - qW(2:end)*(obj.Rwall{1+Np, 2} + ...
                                  obj.Rcw{1+Np, 2})/(obj.T0 - obj.Tinf);
            qWall_(:, 1) = qW/(2*pi*obj.rWL(Np+1));
            for i = 2:N
                thetaW_(2:end, i) = thetaW_(2:end, i-1) - ...
                                    qW(2:end)*(obj.Rwall{i+Np, 2} + ...
                                    obj.Rcw{i+Np, 2})/(obj.T0 - obj.Tinf);
                qWall_(:, i) = qW/(2*pi*obj.rWL(Np+i));
            end 
        end
        function setLumpedParticleWallMesh(obj)
            % generates a mesh for the particle domain and the specified
            % wall insulation layers
            N = size(obj.wallInsulation, 1);
            obj.meshedWallInsulation = {};
            % set particle mesh
            r_ = linspace(0, obj.b*obj.Hp, obj.nrbar+1);
            for j = 1:obj.nrbar
                obj.meshedWallInsulation{j, 1} = 'particles';
                obj.meshedWallInsulation{j, 2} = [r_(j), r_(j+1)];
                obj.meshedWallInsulation{j, 3} = obj.kp;
                obj.meshedWallInsulation{j, 4} = obj.rhopPack;
                obj.meshedWallInsulation{j, 5} = obj.cpp;
                obj.meshedWallInsulation{j, 6} = 1e16;
                obj.thetaWall(j) = 1;                                
            end
            p = obj.nrbar+1;
            for i = 1:N
                M = obj.wallInsulation{i, 6};
                rmin = obj.wallInsulation{i, 2}(1);
                rmax = obj.wallInsulation{i, 2}(2);
                r_ = linspace(rmin, rmax, M+1);
                for j = 1:M
                    obj.meshedWallInsulation{p, 1} = ...
                        obj.wallInsulation{i, 1};
                    obj.meshedWallInsulation{p, 2} = [r_(j), r_(j+1)];
                    obj.meshedWallInsulation{p, 3} = ...
                        obj.wallInsulation{i, 3};
                    obj.meshedWallInsulation{p, 4} = ...
                        obj.wallInsulation{i, 4};
                    obj.meshedWallInsulation{p, 5} = ...
                        obj.wallInsulation{i, 5};
                    if j == 1
                        obj.meshedWallInsulation{p, 6} = ...
                        obj.wallInsulation{i, 7};
                    else
                        obj.meshedWallInsulation{p, 6} = 1e16;
                    end
                    obj.thetaWall(p) = obj.wallInsulation{i, 8};
                    p = p+1;
                end             
            end 
            obj.rWL = unique(cell2mat(obj.meshedWallInsulation(:, 2)), 'sorted');
            obj.rWL = obj.rWL(1:end-2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % internal air and roof heat transfer functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initializeER(obj)
            % initialized state variables for enclosure radiation system
            N = size(obj.roofInsulation, 1); M = size(obj.wallInsulation, 1);
            setERMassMatrix(obj);
            setERTimeConstants(obj);
            obj.thetaRoof = zeros(1, N);
            obj.thetaWA = zeros(1, M);
            y0 = computeERInitialCondition(obj);
            obj.qRadR = qs2QR(obj, y0(N+M+2));
            obj.qRadW = qs2QW(obj, y0(N+M+3));            
        end        
        function computeER(obj, Fo_)
            % computes the temperatures and heat transfer data for the
            % enclosure with the current model data
            N = size(obj.roofInsulation, 1); M = size(obj.wallInsulation, 1);
            setERTimeConstants(obj);
            y0 = computeERInitialCondition(obj);
            options = odeset('Mass', obj.Mass, 'RelTol', 1e-6);
            [~, y] = ode15s(@(t, y) erDAE(obj, t, y), Fo_, y0, options);
            obj.thetaRoof = y(end, 1:N);
            obj.thetaWA = y(end, N+1:N+M);
            obj.thetaA = y(end, N+M+1);
            obj.qRadR = qs2QR(obj, y(end, N+M+2));
            obj.qRadW = qs2QW(obj, y(end, N+M+3));         
        end
        function setERTimeConstants(obj)
            % computes the time constants for the roof and wall odes
            N = size(obj.roofInsulation, 1); M = size(obj.wallInsulation, 1);
            obj.Rroof = {}; obj.CapRoof = {}; obj.ResWA = {}; obj.CapWA = {};
            % roof insulation layer capacitance and resistance
            for i = 1:N
                obj.Rroof{i, 1} = obj.roofInsulation{i, 1};
                obj.CapRoof{i, 1} = obj.roofInsulation{i, 1};
                t_ = abs((obj.roofInsulation{i, 2}(2) - ...
                    obj.roofInsulation{i, 2}(1)));
                k_ = obj.roofInsulation{i, 3};
                rho_ = obj.roofInsulation{i, 4};
                c_ = obj.roofInsulation{i, 5};
                obj.Rroof{i, 2} = t_/(pi*(obj.bp*obj.Hp)^2*k_);
                obj.CapRoof{i, 2} = rho_*c_*t_*pi*(obj.bp*obj.Hp)^2;        
            end    
            % convective resistance
            obj.Rroof{N + 1, 1} = 'external convection';
            obj.Rroof{N + 1, 2} = 1/(pi*(obj.bp*obj.Hp)^2*obj.hInf); 
            obj.Rroof{N + 2, 1} = 'internal convection';
            obj.Rroof{N + 2, 2} = 1/(pi*(obj.bp*obj.Hp)^2*obj.ha);
            obj.tauRA = obj.Hp^2/(obj.alphapPacked*obj.CapRoof{1, 2}*obj.Rroof{N+2, 2});
            % time constants
            obj.tauR = NaN*ones(N, 2);
            obj.tauR(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                              [obj.CapRoof{:, 2}].*[obj.Rroof{1:end-2, 2}]);
            obj.tauR(:, 2) = [obj.Hp^2./(obj.alphapPacked* ...
                                [obj.CapRoof{2:end, 2}].*[obj.Rroof{1:end-3, 2}]), NaN];
            obj.tauRinf = obj.Hp^2./(obj.alphapPacked* ...
                                obj.CapRoof{end, 2}*(obj.Rroof{N, 2} + obj.Rroof{N+1, 2}));
            % wall insulation layer capacitance and resistance
            for i = 1:M
                obj.ResWA{i, 1} = obj.wallInsulation{i, 1};
                obj.CapWA{i, 1} = obj.wallInsulation{i, 1};
                r1 = obj.wallInsulation{i, 2}(1);
                r2 = obj.wallInsulation{i, 2}(2);
                k_ = obj.wallInsulation{i, 3};
                rho_ = obj.wallInsulation{i, 4};
                c_ = obj.wallInsulation{i, 5};
                obj.ResWA{i, 2} = log(r2/r1)/(2*pi*obj.Hp*(1-obj.ztop)*k_);
                obj.CapWA{i, 2} = rho_*c_*obj.Hp*(1-obj.ztop)*pi*(r2^2 - r1^2);
            end
            % convective resistance
            obj.ResWA{M + 1, 1} = 'external convection';
            obj.ResWA{M + 1, 2} = 1/(2*pi*r2*obj.Hp*(1-obj.ztop)*obj.hInf);
            obj.ResWA{M + 2, 1} = 'internal convection';
            obj.ResWA{M + 2, 2} = 1/(2*pi*obj.bp*obj.Hp*obj.Hp*(1-obj.ztop)*obj.ha);
            obj.tauWA = obj.Hp^2/(obj.alphapPacked*obj.CapWA{1, 2}*obj.ResWA{M+2, 2});
            % time constants
            obj.tauW = NaN*ones(M, 2);
            obj.tauW(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                              [obj.CapWA{:, 2}].*[obj.ResWA{1:end-2, 2}]);
            obj.tauW(:, 2) = [obj.Hp^2./(obj.alphapPacked* ...
                                [obj.CapWA{2:end, 2}].*[obj.ResWA{1:end-3, 2}]), NaN];
            obj.tauWinf = obj.Hp^2./(obj.alphapPacked* ...
                                obj.CapWA{end, 2}*(obj.ResWA{M, 2} + obj.ResWA{M+1, 2}));
            % air temperature model
            computeAirSpecificHeat(obj);
            obj.ResP = {}; obj.CapAir = {};
            obj.ResP{1} = 'internal convection';
            obj.ResP{2} = 1/(pi*(obj.bp*obj.Hp)^2*obj.ha);
            obj.CapAir{1} = 'internal air';
            obj.CapAir{2} = obj.rhopPack*obj.cA*pi*(obj.bp*obj.Hp)^2*obj.Hp*(1 - obj.ztop);
            obj.tauAR = obj.Hp^2/(obj.alphapPacked*obj.CapAir{2}*obj.Rroof{N+2, 2});
            obj.tauAW = obj.Hp^2/(obj.alphapPacked*obj.CapAir{2}*obj.ResWA{M+2, 2});
            obj.tauAP = obj.Hp^2/(obj.alphapPacked*obj.CapAir{2}*obj.ResP{2});                                    
        end 
        function setERMassMatrix(obj)
            N = size(obj.roofInsulation, 1); M = size(obj.wallInsulation, 1);
            obj.Mass = eye(N+M+1);
            obj.Mass = [obj.Mass, zeros(N+M+1, 2)];
            obj.Mass = [obj.Mass; zeros(2, N+M+3)];
            obj.Mass = sparse(obj.Mass);
        end
        function y0 = computeERInitialCondition(obj)
            % computes the initial condition for all state variables
            yR1_ = obj.thetaRoof(1);
            yW1_ = obj.thetaWA(1);
            y0_ = erInitialConditionFcn(obj, yR1_, yW1_); 
            y0 = [obj.thetaRoof'; obj.thetaWA'; obj.thetaA; y0_];
        end
        function y = erInitialConditionFcn(obj, yR1, yW1)
            % computes the initial condition for all states with the
            % initial temperature condition for the wall and roof
            y(1, 1) = obj.H^2/(obj.alphaPacked*obj.CapRoof{1, 2}*(obj.H/obj.Hp)^3*(obj.T0 - obj.Tinf)) ...
                    *(((obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Frw*obj.epsilonR ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.Frw*obj.epsilonR ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonR ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP ...
                    + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Frw*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.Frw*obj.epsilonW ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    + obj.Ap*obj.Aw*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.epsilonP*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Aw*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW) ...
                    ./ ...
                    (obj.Ar^2*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    + obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP ...
                    + obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonR ...
                    + obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonR ...
                    + obj.Ar*obj.Aw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar^2*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonR ...
                    + obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP ...
                    + obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP ...
                    + obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    + obj.Ap*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpr*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW)) ...
                    - obj.sigma*(obj.theta2TK(yR1)).^4) ...
                    ./((1 - obj.epsilonR)/(obj.epsilonR*obj.Ar));
           y(2, 1) = obj.H^2/(obj.alphaPacked*obj.CapWA{1, 2}*(obj.H/obj.Hp)^3*(obj.T0 - obj.Tinf)) ...
                    *(((obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Frw*obj.epsilonR ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.Frw*obj.epsilonR ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonR ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP ...
                    + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Frw*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.Frw*obj.epsilonW ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.epsilonP*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar^2*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(yR1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(yW1)).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW) ...
                    ./ ...
                    (obj.Ar^2*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    + obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP ...
                    + obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonR ...
                    + obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonR ...
                    + obj.Ar*obj.Aw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar^2*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonR ...
                    + obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP ...
                    + obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP ...
                    + obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    + obj.Ap*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpr*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                    - obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Ar*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ap*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    - obj.Ar*obj.Aw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                    + obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW)) ...
                    - obj.sigma*(obj.theta2TK(yW1)).^4) ...
                    ./((1 - obj.epsilonW)/(obj.epsilonW*obj.Aw));                         
        end
        function dae = erDAE(obj, ~, y)
            % constructs the system of daes used to solve for the state
            % variables in the radiation enclosure system
            N = size(obj.roofInsulation, 1); M = size(obj.wallInsulation, 1);
            dae = NaN*ones(N+M+3, 1);
            for i = 1:N+M+3
               if i == 1
                    dae(i) = y(N+M+2, :) - obj.tauR(i, 1)*(y(i, :) - ...
                        y(i+1, :)) + obj.tauRA*(y(N+M+1, :) - y(i, :));
               elseif i > 1 && i < N
                    dae(i) = obj.tauR(i-1, 2)*(y(i-1, :) - y(i, :)) - ...
                             obj.tauR(i, 1)*(y(i, :) - y(i+1, :));
               elseif i == N
                    dae(i) = obj.tauR(i-1, 2)*(y(i-1, :) - y(i, :)) - ...
                             obj.tauRinf*y(i, :);
               elseif i == N+1
                    dae(i) = y(N+M+3, :) - obj.tauW(i-N, 1)*(y(i, :) - ...
                        y(i+1, :)) + obj.tauWA*(y(N+M+1, :) - y(i, :));
               elseif i > N+1 && i < N+M
                    dae(i) = obj.tauW(i-N-1, 2)*(y(i-1, :) - y(i, :)) - ...
                             obj.tauW(i-N, 1)*(y(i, :) - y(i+1, :));
               elseif i == N+M
                    dae(i) = obj.tauW(i-N-1, 2)*(y(i-1, :) - y(i, :)) - ...
                             obj.tauWinf*y(i, :);
               elseif i == N+M+1
                    dae(i) = obj.tauAP*(obj.T2thetaK(obj.Tp) - y(i, :)) - ...
                        obj.tauAR*(y(i, :) - y(1, :)) - obj.tauAW*(y(i, :) - y(N+1, :));
               elseif i == N+M+2
                    dae(i) = y(i, :) ...
                        - obj.H^2/(obj.alphaPacked*obj.CapRoof{1, 2}*(obj.H/obj.Hp)^3*(obj.T0 - obj.Tinf)) ...
                        *(((obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonR ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonR ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonR ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP ...
                        + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonW ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        + obj.Ap*obj.Aw*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.epsilonP*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Aw*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW) ...
                        ./ ...
                        (obj.Ar^2*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        + obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP ...
                        + obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonR ...
                        + obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonR ...
                        + obj.Ar*obj.Aw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar^2*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonR ...
                        + obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP ...
                        + obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP ...
                        + obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        + obj.Ap*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpr*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW)) ...
                        - obj.sigma*(obj.theta2TK(y(1, :))).^4) ...
                        ./((1 - obj.epsilonR)/(obj.epsilonR*obj.Ar));
               elseif i == N+M+3
                    dae(i) = y(i, :) ...
                        - obj.H^2/(obj.alphaPacked*obj.CapWA{1, 2}*(obj.H/obj.Hp)^3*(obj.T0 - obj.Tinf)) ...
                        *(((obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonR ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonR ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonR ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP ...
                        + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonW ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.epsilonP*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap^2*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar^2*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.theta2TK(y(1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*(obj.sigma*(obj.Tp)^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*(obj.sigma*(obj.theta2TK(y(N+1, :))).^4)*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW) ...
                        ./ ...
                        (obj.Ar^2*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        + obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP ...
                        + obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonR ...
                        + obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonR ...
                        + obj.Ar*obj.Aw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar^2*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonR ...
                        + obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP ...
                        + obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP ...
                        + obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        + obj.Ap*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpr*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR ...
                        - obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Ar*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ap*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpr*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        - obj.Ar*obj.Aw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap^2*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar^2*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar^2*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Ar*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ap*obj.Aw*obj.Fpr*obj.Fpw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpr*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW ...
                        + obj.Ar*obj.Aw*obj.Fpw*obj.Frw*obj.epsilonP*obj.epsilonR*obj.epsilonW)) ...
                        - obj.sigma*(obj.theta2TK(y(N+1, :))).^4) ...
                        ./((1 - obj.epsilonW)/(obj.epsilonW*obj.Aw));              
               end                
            end            
        end
        function computeAirSpecificHeat(obj)
            % computes the specific heat of air as a function of the
            % internal air temperature (in degC)
            T_ = [-23.15; 26.85; 76.85; 126.85; 176.85; 226.85; 276.85; ...
                326.85; 376.85; 426.85; 476.85; 526.85; 626.85; 726.85; ...
                826.85; 926.85; 1026.85; 1126.85; 1226.85]; % degC
            cp_ = [1.003; 1.005; 1.008; 1.013; 1.020; 1.029; 1.040; ...
                1.051; 1.063; 1.075; 1.087; 1.099; 1.121; 1.142; 1.155; ...
                1.173; 1.19; 1.204; 1.216]*10^3; %J/kgK
            obj.cA = interp1(T_, cp_, obj.theta2T(obj.thetaA));          
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1D spatially resolved transient wall and particle model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeCWP1D(obj, Fo_)
            % computes the continuous 1D temperature solution in the wall
            % with the particle domain either for a single timestep or for
            % an entire time horizon (possible if the initial condition is
            % trivial)
            % initialize wall mesh and initial condition
            buildWallCWP1D(obj);
            % compute eigenvalues and coefficients
            computeCWm(obj, 0);
            % compute temperature distribution for all times in Fo_
            M = size(obj.wallInsulation, 1);            
            for i = 1:length(obj.CWm)                
                for n = 1:M+1
                    bWm_ = obj.bWm{i};
                    if n == 1
                        A = 1; B = 0; alpha_ = obj.alphapPacked;
                    else
                        A = bWm_(2*n-2); B = bWm_(2*n-1);
                        alpha_ = [obj.wallInsulation{n-1, 3}]./ ...
                                ([obj.wallInsulation{n-1, 4}].* ...
                                [obj.wallInsulation{n-1, 5}]);
                    end
                    [x, y, ~, ~] = phiW(obj, obj.rbarW{n}, obj.etaW(i), alpha_);
                    Rn = A*x + B*y;
                    F = XWt(obj, Fo_, obj.etaW(i));
                    obj.thetaW1D{n} = obj.thetaW1D{n} + obj.CWm(i)*F*Rn;                                                       
                end                                
            end            
        end
        function f = XWt(~, Fo_, eta_)
            % transient component of analytic solutions
            f = exp(-eta_^2*Fo_);
        end
        function [x, y, xp, yp] = phiW(obj, r_, eta_, alphai)
            % r-dimension eigenfunctions with eigenvalue, eta_, and thermal
            % diffusivity alpha_
            x = besselj(0, sqrt(obj.alphapPacked/alphai)*eta_*r_); 
            y = bessely(0, sqrt(obj.alphapPacked/alphai)*eta_*r_);
            xp = -sqrt(obj.alphapPacked/alphai)*eta_.* ...
                besselj(1, sqrt(obj.alphapPacked/alphai)*eta_*r_); 
            yp = -sqrt(obj.alphapPacked/alphai)*eta_.* ...
                bessely(1, sqrt(obj.alphapPacked/alphai)*eta_*r_);
        end 
        function ni = NWm(obj, eta_, bWm_)
            % r-BVP eigenvalue problem norm
            M = size(obj.wallInsulation, 1);
            ri = unique([obj.wallInsulation{:, 2}])/obj.Hp;
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            ni = 0.5*obj.bp^2*(besselj(0, eta_*obj.bp)^2 + ...
                                besselj(1, eta_*obj.bp)^2); 
            for i = 2:M
                A = bWm_(2*i-2); B = bWm_(2*i-1);
                r1 = ri(i-1); r2 = ri(i); 
                q1 = sqrt(obj.alphapPacked/alphai(i-1))*eta_*r1; 
                q2 = q1*r2/r1;
                ni = ni + 0.5*(r2^2*(A^2*besselj(0, q2)^2 + ...
                        2*A*B*besselj(0, q2)*bessely(0, q2) + ...
                        (A*besselj(1, q2) + B*bessely(1, q2))^2 + ...
                        B^2*bessely(0, q2)^2) - ...
                        r1^2*(A^2*besselj(0, q1)^2 + ...
                        2*A*B*besselj(0, q1)*bessely(0, q1) + ...
                        (A*besselj(1, q1) + B*bessely(1, q1))^2 + ...
                        B^2*bessely(0, q1)^2));
            end
        end 
        function zi = ZWm(obj, eta_)
            % transcendental equation for finding r-dimension eigenvalues
            [A, ~] = computeAWm(obj, eta_); zi = det(A);
            if ~isreal(zi), zi = 1; end         
        end
        function [A, c] = computeAWm(obj, eta_)
            % computes the r-BVP boundary condition matrix
            M = size(obj.wallInsulation, 1);
            A = zeros(2*M + 1);
            c = zeros(2*M, 1);
            ri = unique([obj.wallInsulation{:, 2}])/obj.Hp;
            ki = [obj.wallInsulation{:, 3}];
            rhoi = [obj.wallInsulation{:, 4}];
            cpi = [obj.wallInsulation{:, 5}];
            alphai = ki./(rhoi.*cpi);
            % fill info for particle domain
            Bic = obj.hcw*obj.Hp/obj.kp;
            [x1, ~, xp1, ~] = phiW(obj, ri(1), eta_, obj.alphapPacked);
            [x2, y2, xp2, yp2] = phiW(obj, ri(1), eta_, alphai(1));
            A(1:2, 1:3) = [xp1 + Bic*x1, -Bic*x2, -Bic*y2; obj.kp*xp1, -ki(1)*xp2, -ki(1)*yp2];
            c(1) = -(xp1 + Bic*x1); c(2) = -obj.kp*xp1;
            % fill boundaries at composite connections
            for i = 1:M-1
                Bic = obj.hcw*obj.Hp/ki(i);
                [x1, y1, xp1, yp1] = phiW(obj, ri(i+1), eta_, alphai(i));
                [x2, y2, xp2, yp2] = phiW(obj, ri(i+1), eta_, alphai(i+1));
                A(2*i+1:2*(i+1), 2*i:2*i+3) = [xp1 + Bic*x1, yp1 + Bic*y1, -Bic*x2, -Bic*y2; ...
                     ki(i)*xp1, ki(i)*yp1, -ki(i+1)*xp2, -ki(i+1)*yp2];
            end
            % fill exterior surface boundary
            [x, y, xp, yp] = phiW(obj, ri(end), eta_, alphai(M));
            A(end, end-1:end) = [xp + obj.Biw2*x, yp + obj.Biw2*y];                        
        end
        function c = FourierCoefficientCWm(obj, eta_, bWm_)
            Cnum = mean(mean(obj.rhoW1D{1}))*obj.b*besselj(1, eta_*obj.b)/eta_;
            c = Cnum/NWm(obj, eta_, bWm_);                       
        end
        function computeEtaW(obj, ShowPlot)
            % computes r-BVP eigenvalues, corresponding boundary
            % condition matrices, and eigenfunction coefficients
            interval = linspace(1e-7, obj.qfW, obj.qW);   % interval/spacing 
                                                     % of root calculation
            rm = NaN*ones(obj.qW, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.qW
                rm(i) = fzero(@(eta_) ZWm(obj, eta_), interval(i), options);
            end
            etaW_ = rm(diff(rm) > 1e-10);         % only keep unique roots
            etaW_ = etaW_(etaW_ > 1e-10);   % remove zeros  
            % save eta and eta dependencies for future reference
            obj.etaWStatic = etaW_;
            obj.etaW = etaW_;
            obj.bWm = cell(length(obj.etaW), 1);
%             obj.etaWN = obj.etaW;
            for i = 1:length(obj.etaW)                            
                [A_, c] = computeAWm(obj, obj.etaW(i));
                b_ = A_(1:end-1, 2:end)\c;
                obj.bWm{i} = [1; b_]; 
%                 obj.etaWN(i) = NWm(obj, obj.etaW(i), obj.bWm{i});
            end
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                for i = 1:length(interval)
                    plot(interval(i), ZWm(obj, interval(i)), '.k');
                    hold on
                end
                for i = 1:length(obj.etaW)
                    plot(obj.etaW(i), ZWm(obj, obj.etaW(i)), '.r');
                end
                xlabel('$\hat{\eta}$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$Z_m(\hat{\eta})$', 'interpreter', 'latex',...
                    'FontSize', 14);
                xlim([0, obj.qfW]);
                hold off 
            end    
        end
        function computeCWm(obj, ShowPlot)
            % computes Fourier series coefficients for 1D continuous wall
            % model
            % compute eta values if not already populated
            if isempty(obj.etaWStatic)
                if nargin > 2 && ShowPlot
                    computeEtaW(obj, true);
                else
                    computeEtaW(obj);
                end     
            end           
            if obj.miW > length(obj.etaW) 
                    miW_ = length(obj.etaW); 
            else
                    miW_ = obj.miW;
            end
            CmTemp = NaN*ones(miW_);
            etaTemp = NaN*ones(miW_); 
            bWmTemp = cell(miW_, 1);
            for i = 1:miW_
                CmTemp(i) = FourierCoefficientCWm(obj, obj.etaW(i), ...
                                                            obj.bWm{i});
                etaTemp(i) = obj.etaW(i);
                bWmTemp{i} = obj.bWm{i};                                                
            end
            obj.cGetW = abs(CmTemp) > obj.climW;
            obj.CWm = CmTemp(obj.cGetW);
            obj.etaW = etaTemp(obj.cGetW);
            obj.bWm = bWmTemp(obj.cGetW);
            cGetIdx = find(obj.cGetW);
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', ...
                        'Position', [0 0 0.4 0.25], 'Visible', 'on');
                Cget = zeros(length(obj.CWm), 1);
                for i = 1:length(obj.CWm)
                    Cget(i) = obj.CWm(i);                        
                end
                plot(cGetIdx, Cget, 'rx'); hold on;
                idx = 1:length(CmTemp);
                plot(idx, CmTemp, '-k')
                legend('Used', 'Computed', 'interpreter', 'latex', ...
                    'FontSize', 14, 'NumColumns', 2);
                xlabel('$m$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$CW_{m}(\eta _m)$', 'interpreter', ...
                    'latex', 'FontSize', 14);
                xlim([0, miW_]);
                hold off                 
            end            
        end
        function buildWallCWP1D(obj)
            % assembles wall parameters from the materials defined in
            % wallInsulation
            [obj.rbarW{1}, ~] = nodeGen(obj, [1e-3, obj.b], obj.nrbar);
            obj.thetaW1D{1} = zeros(1, length(obj.rbarW{1}));
            obj.rhoW1D{1} = ones(1, length(obj.rbarW{1}));
            for i = 2:size(obj.wallInsulation, 1) + 1
                [obj.rbarW{i}, ~] = nodeGen(obj, ...
                    obj.wallInsulation{i-1, 2}./obj.Hp, obj.nrbarW{i-1});
                obj.thetaW1D{i} = ...
                    zeros(1, length(obj.rbarW{i}));
                obj.rhoW1D{i} = zeros(1, length(obj.rbarW{i}));
            end
            obj.Biw1 = obj.hcw*obj.Hp/obj.wallInsulation{1, 3};
            obj.Biw1A = obj.hcwA*obj.Hp/obj.wallInsulation{1, 3};
            obj.Biw2 = obj.hInfp*obj.Hp/obj.wallInsulation{end, 3};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conversions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function T = theta2T(obj, theta)
            T = theta*(obj.T0 - obj.Tinf) + obj.Tinf;
        end
        function T = theta2TK(obj, theta)
            T = theta*(obj.T0 - obj.Tinf) + obj.Tinf + 273.15;
        end
        function theta_ = T2theta(obj, T)
            theta_ = (T - obj.Tinf)/(obj.T0 - obj.Tinf);
        end
        function theta_ = T2thetaK(obj, T)
            theta_ = (T - (obj.Tinf+273.15))/(obj.T0 - obj.Tinf);
        end
        function t = Fo2t(obj, Fo, prototype)
            if nargin < 3, prototype = 0; end
            if prototype 
                t = obj.Hp^2*Fo/obj.alphapPacked;
            else
                t = obj.H^2*Fo/obj.alphaPacked;
            end
        end
        function Fo_ = t2Fo(obj, t, prototype)
            if nargin < 3, prototype = 0; end
            if prototype
                Fo_ = t*obj.alphapPacked/obj.Hp^2;
            else
                Fo_ = t*obj.alphaPacked/obj.H^2;
            end                
        end       
        function Q = qs2QW(obj, qs)
            Q = qs*(obj.Hp^2/(obj.alphapPacked*obj.CapWA{1, 2}*(obj.T0 - obj.Tinf)))^-1;
        end
        function Q = qs2QR(obj, qs)
            Q = qs*(obj.Hp^2/(obj.alphapPacked*obj.CapRoof{1, 2}*(obj.T0 - obj.Tinf)))^-1;
        end
        function qs = Q2qsW(obj, Q)
            qs = Q*obj.Hp^2/(obj.alphapPacked*obj.CapWA{1, 2}*(obj.T0 - obj.Tinf));
        end
        function qs = Q2qsR(obj, Q)
           qs = Q*obj.Hp^2/(obj.alphapPacked*obj.CapRoof{1, 2}*(obj.T0 - obj.Tinf));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % general
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [x, dx] = nodeGen(~, xlim, n)
            % generates a set of n chebyshev nodes spaced between xlim(1)
            % and xlim(2)
%             r_ = (xlim(2) - xlim(1))/2; theta_ = linspace(pi, 0, n);
%             x = xlim(1) + r_*(1 + cos(theta_));
%             dx = (eye(n) - diag(ones(1, n-1), -1))*x'; dx = dx(2:end);
              x = linspace(xlim(1), xlim(2), n);
              dx = x(2) - x(1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotting and vizualization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function data = animateRadialLine(obj, rr, plot_filename)
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'RadialLine.gif'; end
            data = cell(length(obj.Fo), 3);         
            figure('Units', 'normalized', 'color', 'white', ...
                'Position', [0 0 0.5 0.3]);
            rw_ = []; thetaw_ = [];
            for i = 1:length(obj.rbarW)
                rw_ = [rw_, obj.rbarW{i}*obj.Hp]; 
                thetaw_ = [thetaw_, obj.thetaW1D{i}(1, :)];
            end
            l1 = plot(rw_, obj.theta2T(thetaw_), '-k', 'LineWidth', 2);
            xlim([min(rw_), max(rw_)])
            ylim([0, 800])
            timeTitle = title( ...
                   sprintf('$t$ = %1.0f h', 0), 'Interpreter', 'latex', 'FontSize', 14);
            xlabel('$r$ ($m$)', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$T$ ($^\circ C$)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...l
                'TickLabelInterpreter', 'latex', 'FontSize', 12)
            data{1, 1} = 0;
            data{1, 2} = rw_;
            data{1, 3} = thetaw_;
            gif(plot_filename, 'frame', gcf);
            for k_ = 2:length(obj.Fo)
                gif;
                t_ = obj.Fo2t(obj.Fo(k_), 1);
                thetaw_ = [];
                for i = 1:length(obj.rbarW)
                    thetaw_ = [thetaw_, obj.thetaW1D{i}(k_, :)];
                end
                timeTitle.String = sprintf('$t$ = %1.0f h', t_/3600);
                set(l1, 'YData', obj.theta2T(thetaw_));
                data{k_, 1} = t_;
                data{k_, 2} = rw_;
                data{k_, 3} = thetaw_;
                pause(rr);
            end
        end       
        function data = animateDiscreteRadialLine(obj, rr, plot_filename, ShowPlot)
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'RadialTempLine.gif'; end
            if nargin < 4, ShowPlot = 1; end
            data = cell(length(obj.Fo), 3);
            data{1, 1} = 0;
            data{1, 2} = obj.rWL;
            data{1, 3} = obj.thetaWall(1, :);
            if ShowPlot
                figure('Units', 'normalized', 'color', 'white', ...
                    'Position', [0 0 0.5 0.3]);            
                l1 = plot(obj.rWL, obj.theta2T(obj.thetaWall(1, :)), '-k', 'LineWidth', 2);
                xlim([min(obj.rWL), max(obj.rWL)])
                ylim([0, 800])
                timeTitle = title( ...
                       sprintf('$t$ = %1.1f h', 0), 'Interpreter', 'latex', 'FontSize', 14);
                xlabel('$r$ ($m$)', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$T$ ($^\circ C$)', 'interpreter', 'latex', 'FontSize', 14);
                set(gca, 'box', 'off', 'TickDir', 'both', ...
                    'TickLength', [0.01, 0.025], ...l
                    'TickLabelInterpreter', 'latex', 'FontSize', 12)
                gif(plot_filename, 'frame', gcf);
            end
            for k_ = 2:length(obj.Fo)
                t_ = obj.Fo2t(obj.Fo(k_), 1); 
                if ShowPlot
                    gif;
                    timeTitle.String = sprintf('$t$ = %1.1f h', t_/3600);
                    set(l1, 'YData', obj.theta2T(obj.thetaWall(k_, :)));
                    pause(rr);
                end               
                data{k_, 1} = t_;
                data{k_, 2} = obj.rWL;
                data{k_, 3} = obj.thetaWall(k_, :);                
            end
        end 
        function [eT, eq] = plotSteadyStateModelError(~, T, Tss, q, qss, t, r, ShowPlot)
            figure('Units', 'normalized', 'color', 'white', ...
                        'Position', [0 0 0.5 0.4], 'visible', ShowPlot);
            eT = 100*abs(trapz(r, T, 2) - trapz(r, Tss, 2))./trapz(r, T, 2);
            eq = 100*abs(trapz(r, q, 2) - trapz(r, qss, 2))./trapz(r, q, 2);
            plot(t/3600, eT, 'LineWidth', 1); hold on;
            plot(t/3600, eq, 'LineWidth', 1);
            legend('Temperature Error', 'Heat Loss Error', ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Percent Error', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14);
            xlim([0, max(t)/3600]);
        end
        function data = compareRadialTemps(obj, data1, data2, rr, plot_filename, data1name, data2name)
            if nargin < 4, rr = 0.5; end
            if nargin < 5, plot_filename = 'RadialTempLine.gif'; end
            if nargin < 6, data1name = 'test1'; end
            if nargin < 7, data2name = 'test2'; end
            data = cell(length(obj.Fo), 3);
            figure('Units', 'normalized', 'color', 'white', ...
                'Position', [0 0 0.5 0.3]);
            Hp_ = 7;
            Hpf_ = 7*14.29;
            l1 = plot(data1{1, 2}./Hp_, obj.theta2T(data1{1, 3}), '-k', 'LineWidth', 2);
            hold on;
            l2 = plot(data2{1, 2}./Hpf_, obj.theta2T(data2{1, 3}), '-r', 'LineWidth', 2);
            xline(obj.b, '--k', 'LineWidth', 1);            
            xlim([min(data1{1, 2}./Hp_), max(data1{1, 2}./Hp_)])
            ylim([0, 800])
            legend(data1name, data2name, 'Particle-Wall Boundary', 'interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
            timeTitle = title( ...
                   sprintf('$t$ = %1.1f h', 0), 'Interpreter', 'latex', 'FontSize', 14);
            xlabel('$\overline{r}$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$T$ ($^\circ C$)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...l
                'TickLabelInterpreter', 'latex', 'FontSize', 12)
            gif(plot_filename, 'frame', gcf);
            for k_ = 2:size(data1, 1)
                gif;
                t_ = data1{k_, 1};
                timeTitle.String = sprintf('$t$ = %1.1f h', t_/3600);
                set(l1, 'YData', obj.theta2T(data1{k_, 3}));
                % set(l1, 'XData', data1{k_, 2});
                hold on;
                % set(l2, 'XData', data2{k_, 2});
                set(l2, 'YData', obj.theta2T(data2{k_, 3}));
                hold on;
                pause(rr);
            end
        end
    end
end

