 classdef FF < handle
    properties   
        % dependent variables
        FoEnd = 1   % simulation end time
        rtop        % r-coordinates for calculations in top boundary
        nrtop = 100 % number of r-nodes in top boundary
        rbar        % r-coordinates for calculations in stagnant region and 
                    % top boundary
        nrbar = 1000 % number of r-nodes to use for stagnant region
        rhat        % r-coordinate for calculations in center channel
        nrhat = 25  % number of r-nodes in center channel
        rbarH       % r-coordinates for 'H' and 'C' computations
        nrH = 600   % number of r-nodes to use for 'H' and 'C' computations
        zcenter     % z-coordinates for calculations in center channel
        nzc0 = 250  % max number of z-nodes in center channel
        nzc         % number of z-nodes in center channel
        zbar        % z-coordinates for calculations in stagnant region
        zbar0       % initial z-coordinates for stagnant region
        nzbar0 = 1500 % initial number of mesh nodes in stagnant region
        nzbar       % number of z-nodes to use for stagnant region
        zhat        % z-coordinates for calculations in the top boundary
        nzhat = 20  % number of z-nodes in top boundary
        zbarH       % z-coordinates for 'H' and 'C' computations
        nzH0 = 2000
        nzH         % number of z-nodes to use for 'H' and 'C' computations
        mubar       % spherical coordinate system angular dimension
        rho         % spherical coordinate system radial dimension
        rho0        % initial spherical radial coordinates
        z           % z-coordinates for whole domain
        r           % r-coordinates for whole domain
        zmc         % z-mesh for mass flow cone
        Fo          % non-dimensional time (Fourier number) vector
        FoMode      % storage mode at each time (H, D, C)
        FoNow = 0   % current iteration non-dimensional time
        FoModeNow   % current cycle mode (H, D, C)
        FoEmpty     % non-dimensional time when tank is completely empty
        tEmpty      % time when tank is completely empty
        % temperatures in celcius
        T0 = 1061           % initial temperature of particles in bin
        Tinf = 20           % ambient temperature
        Tref = 0            % reference temperature
        % eigenvalues
        beta            % stagnant radial temperature solution
        eta             % stagnant z temperature solution
        etaFD           % eigenvalues for discharge filter 1
        betaFD          % eigenvalues for discharge filter 2
        betaH           % eigenvalues for holding temperature solution
        etaH            % eigenvalues for holding temperature solution
        betaFH           % eigenvalues for holding solution filter 1
        betaC           % eigenvalues for charging solution
        etaC            % eigenvalues for charging solution
        betaFC           % eigenvalues for charging solution filter 1
        % Fourier coefficients
        PsiCnm          % homogeneous temperature solution coefficients
        GCnm            % Green's Function coefficients
        FDCn            % Coefficients for stagnant solution filter 1
        FDCm            % Coefficients for stagnant solution filter 2
        CnmH            % Coefficients for holding solution
        FHhatCn         % Coefficients for holding solution filter 1
        CnmC            % Coefficients for charging solution
        FCCn            % Coefficients for charging solution filter 1
        % solution matrices       
        psiS            % stagnant region homogeneous temperature solution         
        thetaO          % centerline temperature at outlet
        thetaOB         % bulk temperature at outlet
        thetaS          % stagnant region temperature solution
        thetaSDot       % stagnant region temperature time derivative
        FS0             % initial condition matrix for stagnant region
        FD              % discharg homogeneaous filtering function
        thetaIC         % initial condition contribution in Green's Formula
        thetaICDot      % IC contribution derivative
        RIC             % radial initial condition component (static)
        thetaBC1        % stagnant temperature contribution from BC1
        thetaBC1Dot     % "" derivative
        g1              % nonhomogeneous time-dependent BC1
        g2              % nonhomogeneous BC2
        thetaBC3        % stagnant temperature contribution from BC3
        thetaBC3Dot     % "" derivative
        g3              % nonhomogeneous time-dependent BC3
        g4              % nonhomogeneous BC4
        g4f             % s.s. filter portion of BC4
        thetaT          % top boundary temperature solution
        expTop          % static matrix exponential for top solution
        expMC           % static matrix exponential for mass flow cone
        topLoss         % temperature loss accross top moving surface
        FT0             % initial temperature for top boundary
        thetaC          % center boundary temperature solution
        scPe            % sensitivity w.r.t. Pe
        scQc            % sensitivity w.r.t. Qc
        sca             % sensitivity w.r.t. a
        FC0             % initial temperature for center boundary
        thetaChat       % boundary condition at zhat = 1 for center channel
        thetaI = 0.5    % offset boundary temperature for stagnant region
        thetaWI = 0.5   % offset boundary temerature for wall transient
        rhoH            % initial condition for holding solution
        rhoMC           % initial condition for cone temperature        
        FH              % filtering function for holding solution
        FCh             % filtering function for charging solution
        thetaH          % holding temperature solution sequence
        thetaCh         % charging temperature solution sequence
        thetaCi = 1     % charging inlet particle temperature
        thetaCp         % prescribed temperature boundary temp for charging
        thetaMC         % charging mixing depth temperature
        theta           % full temperature solution
        thetaK          % full temperature solution for Kevin's Model Sim
        thetaMatchD     % discharge temperature with matched mesh
        thetaMatchH     % charge and holding temperature with matched mesh
        ubar            % top boundary radial velocity solution
        uzc             % velocity in cone for mass flow 
        wbar            % center boundary vertical velocity solution
        ur              % radial velocity component at every cell in bin
        uz              % z velocity component at every cell in bin
        AT              % state matrix for top boundary temperature
        OmegaT1, OmegaT2, OmegaT3, OmegaT4, OmegaT5
        LambdaT1, LambdaT2, LambdaT3
        GT
        AC              % state matrix for center channel temperature
        AMC             % state matrix for mass flow cone
        dAC_Pe          % state partial w.r.t. Pe
        dAC_Qc          % state partial w.r.t. Qc
        dAC_a           % state partial w.r.t. a
        OmegaC1, OmegaC2, OmegaC3, OmegaC4, OmegaC5
        LambdaC1, LambdaC2
        GC
        % storage bin wall linear RC model
        Rwall           % (K/W) cell array storing resistance of wall layers
        Rcw             % (K/W) contact resistance between wall layers
        Rbase           % (K/W) cell array storing resistance of base layers
        Rcb             % (K/W) contact resistance between base layers
        Rroof           % (K/W) cell array storing resistance of top layers
        Rcr             % (K/W) contact resistance between roof layers
        ResP            % (K/W) convective resistance at particle top-air
        ResWA           % (K/W) cell array with resistance of exposed wall
        CapWall         % (J/K) cell array storing capacitance of wall layers
        CapBase         % (J/K) cell array storing capacitance of base layers
        CapRoof         % (J/K) cell array storing capacitance of top layers
        CapAir          % (J/K) thermal capacitance of internal air
        CapWA           % (J/K) cell array for capacitance of exposed wall
        thetaWall       % state variables for wall
        rWL             % (m) radial mesh for discretized lumped wall
        thetaWA         % state variables for exposed wall
        thetaBase       % state variables for base
        zBL             % (m) vertical mesh for discretized lumped base
        thetaRoof       % state variables for top     
        zRL             % (m) vertical mesh for discretized lumped roof
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
        % radiation model parameters
        epsilonP = 0.9            % emissivity of particle bed
        epsilonR = 0.8            % emissivity of exposed roof
        epsilonW = 0.8            % emissivity of exposed wall
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
        % wall and base 2D continuouse model variables
        thetaW          % composite wall solution/s        
        thetaWIC        % IC contribution for thetaW
        thetaWBC        % BC contribution for thetaW
        gW1             % current inner wall boundary condition
        gW1p            % averaged boundary condition with particles
        gW1a            % averaged boundary condition with air
        KWI             % boundary contribution kernel integrand
        KW              % boundary contribution kernel
        IBCW = []       % stored inner wall boundary conditions, f(z, t)
        rhoW            % composite wall initial condition/s
        zbarW           % z-dimension vectors for composite wall
        rbarW           % r-dimension vectors for composite wall
        nzbarW          % number of nodes used in zbarW
        nrbarW          % number of nodes used in rbarW
        AWm             % r-BVP boundary condition matrix
        bWm             % r eigenfunction coefficients (IC contribution)
        bWBCm           % r eigenfunction coefficients (BC contribution)
        cWr             % lhs of reduced r-BVP boundary condition sysftem
        etaW            % r-BVP eigenvalues reduced for IC contribution
        betaW           % z-BVP eigenvalues reduced for IC contribution
        etaWN           % r-BVP norms
        betaWN          % z-BVP norms
        etaWBC          % r-BVP eigenvalues reduced for BC contribution
        betaWBC         % z-BVP eigenvalues reduced for BC contribution
        hcw = 1         % wall-particle boundary contact coefficient
        hcwA = 1        % wall-tank air convection coefficient
        Biw1            % biot number for inner-most composite layer
        Biw1A           % "" for tank air connection
        Biw2            % biot number for outer-most composite layer        
        betaWStatic     % static betaW values
        etaWStatic      % static etaW values
        AWmStatic       % boundary condition matrix for etaWStatic
        bWmStatic       % r-BVP coefficients for etaWStatic
        GWCnm           % wall coefficients for IC contribution
        GWCm            % 1D wall coefficients for IC contribution
        GWBCCnm         % wall coefficients for BC contribution
        tauW1           % time constant for wall boundary ramp function
        % conservation validation
        continuity      % cell-to-cell evaluation of 2D continuity equation
        energy          % cell-to-cell evaluation of 2D energy equation
        % static r-dimensional variables
        RStatic         % r-dimension eigenfunctions for all eta values
        RFDStatic       % r-dimension discharge filter eigenfunctions 
                        % for all eta values
        RC3Static       % r-dimension eigenfunctions used for BC3
        RC4Static       % r-dimension eigenfunctions used for BC4
        etaStatic       % static eta eigenvalues
        etaFDStatic     % static eta eigenvalues for discharge filter
        etaHStatic      % static eta eigenvalues for holding solution
        RD              % r-dimension eigenfunctions for current eta values
        RFD             % r-dimension eigenfunctions for discharge filter
                        % with current eta values
        RC3             % r-dimension eigenfunctions for boundary 
                        % contribution for current eta values
        RC4             % r-dimension eigenfunctions for boundary 
                        % contribution for current eta values
        % inputs from Modelica model for system integration
        modelicaInputs  % inputs from modelica model (Tin, m_dot_in/out)
        % non-dimensional flow rates (should be equal for correct solution)
        Qc              % center channel
        Qt              % top boundary
        % heat loss at domain boundaries
        qTopP           % (W/m2) heat loss from top of particle region
        qWall           % (W/m2) heat loss from outer-most composite layer
        qLossW          % (W/m2) prescribed heat flux at bin wall
        qLossB          % (W/m2) prescribed heat flux at bin base
        qLossT          % (W/m2) prescribed heat flux at bin top
        qLossWp         % (kW) total wall heat loss computed in particles
        qLossBp         % (kW) total base heat loss computed in particles
        qLossTp         % (kW) total top heat loss computed in particles
        % controlled dimensional volumetric flow rate
        Q = 1e-6        % m3/s
        QCh = 1e-6      % m3/s charging flow rate
        simMass         % actual mass contained in the bin at each time step
        thryMass        % theoretical mass contained at each time step
        % experimental data for comparison and fitting (matched, true)
        zbarExp_ = [], zbarExp = []
        rbarExp_ = [], rbarExp = []
        FoExp_ = [], FoExp = []
        thetaExp   
        thetaOExp          % centerline temperature at outlet
        thetaOBExp         % bulk temperature at outlet
        % prototype dynamic and geometric parameters
        Hp = 7                      % (m) height of prototype bin
        bp = 0.3214                 % (m) inner radius of prototype bin
        hp1 = 10                    % (W/m2K) top flow surface h
        hp2 = 10                    % (W/m2K) bottom overall convection 
                                    % coefficient for prototype bin
        hp3 = 10                    % (W/m2K) center flow channel h
        hp4 = 10                    % (W/m2K) wall overall convection 
                                    % coefficient for prototype bin
        hp5 = 1                     % (W/m2K) free surface h
        hp5D                        % "" for discharge mode
        hp5H                        % "" for holding mode
        hp5C                        % "" for chargeing mode
        hpTop                       % (W/m2K) overall heat transfer
                                    % coefficient for the top of the bin
        kp = 0.7                    % (W/mK) particle packed thermal 
                                    % conductivity
        rhopPack = 2000             % (kg/m3) particle packed bulk 
                                    % density
        rhopLoose = 1810            % (kg/m3) particle loose bulk 
                                    % density
        cpp = 1025.965              % (J/kgK) average particle heat 
                                    % capacity 
        alphapPacked = 3.4114e-7    % (m2/s) particle thermal 
                                    % diffusivity 
        alphapLoose = 3.7695e-7     % (m2/s) particle thermal 
                                    % diffusivity 
        nup = 2.5e-8                % (m2/s) dynamic viscosity of 
                                    % particles moving in air 
        Bip1                 % Biot number at boundary 1
        Bip2                 % "" boundary 2
        Bip3                 % "" boundary 3
        Bip4                 % "" boundary 4
        Bip5                 % "" surface of top boundary open to air
        Bip5D                % "" for discharge mode
        Bip5H                % "" for holding mode
        Bip5C                % "" for chargeing mode
        Gap                  % Galilei number
        Prp                  % Prandtl number
        Rep                  % Reynolds number w.r.t. Uinf
        Pep                  % Peclet number w.r.t. Uinf 
        Uinfp                % prototype outlet velocity
        Qp                   % prototype volumetric flow rate
        QChp                 % prototype charging flow rate
        % variables indicating location for discrete storage points
        zStore, rStore, thetaStore
        zExp, rExp, thetaStoreExp
        % static model parameters with default values
        H = 1                       % (m) height of stagnant region
        ztop                        % zbar element of top surface
        z0C = 0.05                  % starting height for full charge
        deltaM = 0.01               % mixing depth for charging mode
        HNow                        % (m) current height of top surface
        HNowC                       % (m) current height of top surface
        dH                          % (m) change in height per time step
        dHCh                        % "" for charging mode
        h = 0.025                   % nondimensional height of top flow
                                    % boundary
        delta = pi/10               % angle of repose for top surface
        a = 0.01                    % inner nondimensional radius at 
                                    % top of center channel
        a0 = 0.01                   % inner radius at outlet
        aOutlet = 0.01              % radius of aperature where particles 
                                    % flow out of bin
        amc = 0.01                  % outlet apperature for mass flow cone
        bmc                         % z-dependent cone redius
        hCone                       % height of mass flow cone
        phia = 0                    % angle of incline for center channel
        mua                         % cos(phia)
        b = 0.4                     % outer nondimensional radius
        bzc                         % outer radius vector in mass flow cone
        h1 = 10                     % (W/m2-K) heat transfer coefficient at 
                                    % boundary 1
        h2 = 0.07                   % "" boundary 2
        h3 = 5                      % "" boundary 3
        h4 = 0.07                   % "" boundary 4
        h5 = 10                     % "" for free surface
        h5D = 10                    % "" for discharge mode
        h5H = 10                    % "" for holding mode
        h5C = 1                     % "" for chargeing mode
        hInf = 10                   % "" heat transfer coefficient to 
                                    % surroundings
        g = 9.80665                 % (m/s2)
        k = 0.7                     % (W/mK) particle packed thermal 
                                    % conductivity
        rhoPack = 2000              % (kg/m3) particle packed bulk 
                                    % density
        rhoLoose = 1810             % (kg/m3) particle loose bulk 
                                    % density
        cp = 1025.965               % (J/kgK) average particle heat 
                                    % capacity 
        alphaPacked                 % (m2/s) particle thermal 
                                    % diffusivity 
        alphaLoose                  % (m2/s) particle thermal 
                                    % diffusivity 
        mup = 2.5*1.81e-5           % (kg/ms) viscosity of particles 
                                    % moving in air(Bicerano, Douglas 
                                    % and Brune, 1999)
        nu                          % (m2/s) dynamic viscosity of 
                                    % particles moving in air 
        CapS                        % (J/m2K) relative thermal capacitance 
                                    % of stagnant region
        epsilon2                    % ratio of base thermal capacitance to 
                                    % stagnant region thermal capacitance
        epsilon4                    % ratio of wall thermal capacitance to 
                                    % stagnant region thermal capacitance
        thetaA = 0                  % ambient temperature inside
        Uinf = 0.02                 % bulk velocity at outlet  
        baseInsulation              % insulation info for base
        meshedBaseInsulation        % meshed insulation info for base
        wallInsulation              % insulation info for wall
        meshedWallInsulation        % meshed insulation info for wall
        roofInsulation              % insulation info for top of bin
        meshedRoofInsulation        % meshed insulation info for roof
        % nondimensional terms
        Bi1                 % Biot number at boundary 1
        Bi2                 % "" boundary 2
        Bi3                 % "" boundary 3
        Bi4                 % "" boundary 4
        Bi5                 % "" surface of top boundary open to air
        Bi5D                % "" for discharge mode
        Bi5H                % "" for holding mode
        Bi5C                % "" for chargeing mode
        Ga                  % Galilei number
        Pr                  % Prandtl number
        Re                  % Reynolds number w.r.t. Uinf
        Pe                  % Peclet number w.r.t. Uinf  
        % discharge fourier summation parameters
        clim = 5e-7      % only fourier coefficients > clim are used
        cGet             % index array for obtaining eigenvalues
        p = 2000         % total number of beta values computed
        pf = 1000        % range for beta values to be computed
        q = 4000         % total number of eta values computed
        qf = 2500        % range for eta values to be computed
        ni = 100         % number of beta values used in computation of Cnm
        mi = 50          % number of eta values used in computation of Cnm
        climFD = 5e-7    % only fourier coefficients > clim are used
        cGetFD           % index array for obtaining eigenvalues
        pFD = 4000       % total number of beta values computed
        pfFD = 2500      % range for beta values to be computed
        qFD = 4000       % total number of eta values computed
        qfFD = 2500      % range for eta values to be computed
        niFD = 30        % number of beta values used in computation of Cnm
        miFD = 30        % number of eta values used in computation of Cnm
        % holding/mass flow fourier summation parameters
        climH = 5e-7     % only fourier coefficients > clim are used
        cGetH            % index array for obtaining eigenvalues
        pH = 2000        % total number of beta values computed
        pfH = 1000       % range for beta values to be computed
        qH = 4000        % total number of eta values computed
        qfH = 2500       % range for eta values to be computed
        niH = 100        % number of beta values used in computation of Cnm
        miH = 50         % number of eta values used in computation of Cnm
        climFH = 5e-7    % only fourier coefficients > clim are used
        cGetFH           % index array for obtaining eigenvalues
        pFH = 1000       % total number of beta values computed
        pfFH = 100       % range for beta values to be computed
        niFH = 30        % number of beta values used in computation of Cnm
        % wall conduction fourier summation parameters
        climW = 5e-6     % only fourier coefficients > clim are used
        climWBC = 5e-6   % only fourier coefficients > clim are used
        cGetW            % index array for obtaining eigenvalues
        cGetWBC          % index array for obtaining eigenvalues
        pW = 2000        % total number of beta values computed
        pWBC = 2000      % total number of beta values computed
        pfW = 1000       % range for beta values to be computed
        pfWBC = 1000     % range for beta values to be computed
        qW = 2000        % total number of eta values computed
        qWBC = 2000      % total number of eta values computed
        qfW = 1000       % range for eta values to be computed
        qfWBC = 1000     % range for eta values to be computed
        niW = 50         % number of beta values used in computation of Cnm
        niWBC = 50       % number of beta values used in computation of Cnm
        niWBCtot = 10    % number of beta values used in summation
        miW = 50         % number of eta values used in computation of Cnm
        miWBC = 50       % number of eta values used in computation of Cnm
        miWBCtot = 10    % number of eta values used in summation
        % discretization parameters
        df = 0.1        % nondimensional time-step
        dt              % (s) time-step
        dr = 0.005      % radial mesh size for full domain computations
        dz = 0.005      % z-mesh size for full domain computations
        drH = 0.001     % nondimensional radius mesh for 'H'  
        drtop = 0.05    % nondimensional radius top boundary mesh size
        drbar = 0.001   % nondimensional radius large mesh size
        drhat = 0.01    % nondimensional radius small mesh size
        dzH = 0.001     % nondimensional height mesh for 'H'
        dzc = 0.05      % nondimensional height center channel mesh size
        dzbar = 0.001   % nondimensional height large mesh size         
        dzhat = 0.01    % nondimensional height small mesh size
        dzmc = 0.01     % nondimensional mass flow cone mesh size
        modZH           % number of expansions/compressions to perform 
                        % before increasing/decreasing size of z-mesh for
                        % charging/mass-flow discharge mode
        modZS           % number of compressions to perform before 
                        % reducing size of z-mesh for stagnant region
        modZC           % number of compressions to perform before 
                        % reducing size of z-mesh for center region
        dmu = 0.01      % nondimensional spherical angle mesh size
        drho = 0.1      % nondimensional spherical radial mesh size
        % data saving parameters
        ls = 10         % max storage size for time steps
        thetaFolder     % folder to save theta matrices in
        objFolder       % folder to save FF object  
        % data storage for matlab app
        chargeDurration = 0
        holdDurration = 0
        dischargeDurration = 0
        timeStep = 1
        inletMassFlowRate = 0
        outletMassFlowRate = 0
        currentStatus = ''
        % figures and tables
        vtbl            % table containing all static variables
        cfig            % figure showing Fourier coefficients
        betafig         % figure showing beta values
        etafig          % figure showing eta values
        psifig          % figure showing temperature for psi
        thetafig        % figure showing temperature for theta (discharge)
        thetaHfig       % figure showing temperature for theta (holding)
        thetaCfig       % figure showing temperature for theta (charging)
        ubarfig         % figure showing velocity in top boundary
        ubarContFig     % figure showing continuity solution in top bound.
        wbarfig         % figure showing velocity in center boundary
        computeBUfig    % figure showing approximated bulk velocities
        thetaOfig       % figure showing bulk outlet temperature
        betaHfig        % figure showing holding beta eigenvalues
        etaHfig         % figure showing holding eta eigenvalues
        etaFHfig        % figure showing holding filter eta eigenvalues
        cfigFHhat       % figure showing holding filter coefficients
        cfigH           % figure showing holding Fourier coefficients
        betaCfig        % figure showing charging beta eigenvalues
        etaCfig         % figure showing charging eta eigenvalues
        etaFCfig        % figure showing charging filter eta eigenvalues
        cfigFC          % figure showing charging filter coefficients
        cfigC           % figure showing charging Fourier coefficients
    end      
    methods 
        function obj = FF(a_, b_, h_, FoEnd_)
            if nargin > 1
                obj.a = a_;
                obj.b = b_;
                obj.h = h_;
                obj.FoEnd = FoEnd_;
            end 
            obj.ztop = 1;
            [obj.rtop, obj.drtop] = nodeGen(obj, [obj.a, obj.b], obj.nrtop);
            [obj.rbar, obj.drbar] = ... 
                        nodeGen(obj, [obj.a0, obj.b], obj.nrbar);
            [obj.rbarH, obj.drH] = ... 
                        nodeGen(obj, [0, obj.b], obj.nrH);
            [obj.rhat, obj.drhat] = nodeGen(obj, [1e-6, obj.a0], obj.nrhat);
            obj.r = [obj.rhat, obj.rbar(2:end)];
            obj.nzbar = ceil(obj.nzbar0*obj.ztop);
            obj.zbar0 = nodeGen(obj, [0, obj.ztop], obj.nzbar);
            [obj.zbar, obj.dzbar] = ... 
                        nodeGen(obj, [0, obj.ztop], obj.nzbar);
            obj.nzH = ceil(obj.nzH0*obj.ztop);
            [obj.zbarH, obj.dzH] = ... 
                        nodeGen(obj, [0, obj.ztop], obj.nzH);
            obj.nzc = ceil(obj.nzc0*obj.ztop);
            [obj.zcenter, obj.dzc] = nodeGen(obj, [0, obj.ztop], obj.nzc);
            [obj.zhat, obj.dzhat] = nodeGen(obj, [0, obj.h], obj.nzhat);
            obj.z = [obj.zbar, 1+obj.zhat(2:end)];
            obj.zmc = 0:obj.dzmc:obj.hCone;
            obj.HNow = obj.zbar(end)*obj.H;    
            obj.HNowC = obj.zbarH(end)*obj.H;
            % radiation params
            obj.Ap = pi*(obj.bp*obj.Hp)^2;
            obj.Ar = obj.Ap;
            obj.Aw = 2*pi*obj.bp*obj.Hp*(obj.Hp*(1 - obj.ztop));
            r_ = obj.bp/(1 - obj.ztop);
            x = 2 + 1/r_^2;
            y = sqrt(x^2 - 4);
            obj.Fpr = 0.5*(x - y);
            obj.Fpw = 1 - obj.Fpr;
            obj.Frw = obj.Fpw;
            % scaling and property parameters
            obj.alphapPacked = obj.kp/(obj.rhopPack*obj.cpp); 
            obj.alphapLoose = obj.kp/(obj.rhopLoose*obj.cpp);
            obj.alphaPacked = obj.k/(obj.rhoPack*obj.cp); 
            obj.alphaLoose = obj.k/(obj.rhoLoose*obj.cp);
%             obj.CapS = obj.rhoPack*obj.cp*obj.H;
%             obj.epsilon2 = obj.CapBase/obj.CapS;
%             obj.epsilon4 = obj.CapWall/obj.CapS;
            obj.nup = obj.mup/obj.rhoLoose;
            obj.nu = obj.nup;
            obj.Bi1 = obj.h1*obj.H/obj.k;          
            obj.Bi2 = obj.h2*obj.H/obj.k;               
            obj.Bi3 = obj.h3*obj.H/obj.k;               
            obj.Bi4 = obj.h4*obj.H/obj.k;                    
            obj.Bi5D = obj.h5D*obj.H/obj.k;  
            obj.Bi5H = obj.h5H*obj.H/obj.k; 
            obj.Bi5C = obj.h5C*obj.H/obj.k; 
            obj.Ga = obj.g*obj.H^3/obj.nu^2;                
            obj.Pr = obj.nu/obj.alphaPacked;                        
%             obj.Uinf = obj.Q/(pi*(obj.a0*obj.H)^2);
%             obj.Q = obj.Uinf*pi*(obj.a0*obj.H)^2;
%             obj.Re = obj.H*obj.Uinf/obj.nu;           
%             obj.Pe = obj.Re*obj.Pr;
%             obj.Qc = obj.Q/(obj.H^2*obj.Uinf);
            obj.Bip1 = obj.hp1*obj.Hp/obj.kp;          
            obj.Bip2 = obj.hp2*obj.Hp/obj.kp;               
            obj.Bip3 = obj.hp3*obj.Hp/obj.kp;               
            obj.Bip4 = obj.hp4*obj.Hp/obj.kp;                    
            obj.Bip5D = obj.hp5D*obj.Hp/obj.kp;  
            obj.Bip5H = obj.hp5H*obj.Hp/obj.kp; 
            obj.Bip5C = obj.hp5C*obj.Hp/obj.kp;   
            obj.Gap = obj.g*obj.Hp^3/obj.nup^2;                
            obj.Prp = obj.nup/obj.alphapPacked;                        
            obj.Uinfp = obj.Qp/(pi*(obj.a0*obj.Hp)^2);
            obj.Rep = obj.Hp*obj.Uinfp/obj.nup;           
            obj.Pep = obj.Rep*obj.Prp;
            obj.dt = obj.Fo2t(obj.df);
            obj.dH = obj.Q*obj.dt/(pi*(obj.b*obj.H)^2);
            obj.dHCh = obj.QCh*obj.dt/(pi*(obj.b*obj.H)^2);
            obj.modZH = min(ceil(obj.dzH./(obj.dH/obj.H)));
            obj.modZS = min(ceil(obj.dzbar./(obj.dH/obj.H)));
            obj.modZC = min(ceil(obj.dzc./(obj.dH/obj.H)));
            obj.tEmpty = obj.H*pi*(obj.b*obj.H)^2/obj.Q;
            obj.FoEmpty = obj.t2Fo(obj.tEmpty);
            obj.Fo = 0:obj.df:obj.FoEnd;
            obj.FoMode = cell(1, length(obj.Fo));
            % initialize to holding
            obj.FoMode(:) = {'H'};
        end
        function reInitObj(obj)
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.  
            [obj.rtop, obj.drtop] = nodeGen(obj, [obj.a, obj.b], obj.nrtop);
            [obj.rbar, obj.drbar] = ... 
                        nodeGen(obj, [obj.a0, obj.b], obj.nrbar);
            [obj.rbarH, obj.drH] = ... 
                        nodeGen(obj, [0, obj.b], obj.nrH);
            [obj.rhat, obj.drhat] = nodeGen(obj, [1e-6, obj.a0], obj.nrhat);
            obj.r = [obj.rhat, obj.rbar(2:end)];
            obj.nzbar = ceil(obj.nzbar0*obj.ztop);
            obj.zbar0 = nodeGen(obj, [0, obj.ztop], obj.nzbar);
            [obj.zbar, obj.dzbar] = ... 
                        nodeGen(obj, [0, obj.ztop], obj.nzbar);
            obj.nzH = ceil(obj.nzH0*obj.ztop);
            [obj.zbarH, obj.dzH] = ... 
                        nodeGen(obj, [0, obj.ztop], obj.nzH);
            obj.nzc = ceil(obj.nzc0*obj.ztop);
            [obj.zcenter, obj.dzc] = nodeGen(obj, [0, obj.ztop], obj.nzc);
            [obj.zhat, obj.dzhat] = nodeGen(obj, [0, obj.h], obj.nzhat);
            obj.z = [obj.zbar, 1+obj.zhat(2:end)];
            obj.zmc = 0:obj.dzmc:obj.hCone;
            obj.HNow = obj.zbar(end)*obj.H;    
            obj.HNowC = obj.zbarH(end)*obj.H;
            % radiation params
            obj.Ap = pi*(obj.bp*obj.Hp)^2;
            obj.Ar = obj.Ap;
            obj.Aw = 2*pi*obj.bp*obj.Hp*(obj.Hp*(1 - obj.ztop));
            r_ = obj.bp/(1 - obj.ztop);
            x = 2 + 1/r_^2;
            y = sqrt(x^2 - 4);
            obj.Fpr = 0.5*(x - y);
            obj.Fpw = 1 - obj.Fpr;
            obj.Frw = obj.Fpw;
            % scaling and property parameters
            obj.alphapPacked = obj.kp/(obj.rhopPack*obj.cpp); 
            obj.alphapLoose = obj.kp/(obj.rhopLoose*obj.cpp); 
            obj.alphaPacked = obj.k/(obj.rhoPack*obj.cp); 
            obj.alphaLoose = obj.k/(obj.rhoLoose*obj.cp);  
            obj.nup = obj.mup/obj.rhoLoose;
%             obj.CapS = obj.rhoPack*obj.cp*obj.H;
%             obj.epsilon2 = obj.CapBase/obj.CapS;
%             obj.epsilon4 = obj.CapWall/obj.CapS;
            obj.Bi1 = obj.h1*obj.H/obj.k;          
            obj.Bi2 = obj.h2*obj.H/obj.k;               
            obj.Bi3 = obj.h3*obj.H/obj.k;               
            obj.Bi4 = obj.h4*obj.H/obj.k;                    
            obj.Bi5D = obj.h5D*obj.H/obj.k;  
            obj.Bi5H = obj.h5H*obj.H/obj.k; 
            obj.Bi5C = obj.h5C*obj.H/obj.k; 
            obj.Ga = obj.g*obj.H^3/obj.nu^2;                
            obj.Pr = obj.nu/obj.alphaPacked;   
%             obj.Uinf = obj.Q/(pi*(obj.a0*obj.H)^2);
%             obj.Q = obj.Uinf*pi*(obj.a0*obj.H)^2;
            obj.Qc = obj.Q/(obj.H^2*obj.Uinf);
            obj.Re = obj.H*obj.Uinf/obj.nu;                 
            obj.Pe = obj.Re*obj.Pr;
            obj.Bip1 = obj.hp1*obj.Hp/obj.kp;          
            obj.Bip2 = obj.hp2*obj.Hp/obj.kp;               
            obj.Bip3 = obj.hp3*obj.Hp/obj.kp;               
            obj.Bip4 = obj.hp4*obj.Hp/obj.kp;                    
            obj.Bip5D = obj.hp5D*obj.Hp/obj.kp;  
            obj.Bip5H = obj.hp5H*obj.Hp/obj.kp; 
            obj.Bip5C = obj.hp5C*obj.Hp/obj.kp;  
            obj.Gap = obj.g*obj.Hp^3/obj.nup^2;                
            obj.Prp = obj.nup/obj.alphapPacked;                        
            obj.Uinfp = obj.Qp/(pi*(obj.a0*obj.Hp)^2);
            obj.Rep = obj.Hp*obj.Uinfp/obj.nup;           
            obj.Pep = obj.Rep*obj.Prp;
            obj.dt = obj.Fo2t(obj.df);
            obj.dH = obj.Q*obj.dt/(pi*(obj.b*obj.H)^2);
            obj.dHCh = obj.QCh*obj.dt/(pi*(obj.b*obj.H)^2);
            obj.modZH = min(ceil(obj.dzH./(obj.dH/obj.H)));
            obj.modZS = min(ceil(obj.dzbar./(obj.dH/obj.H)));
            obj.modZC = min(ceil(obj.dzc./(obj.dH/obj.H)));
            obj.tEmpty = obj.H*pi*(obj.b*obj.H)^2/obj.Q;
            obj.FoEmpty = obj.t2Fo(obj.tEmpty);
            obj.Fo = 0:obj.df:obj.FoEnd;
        end     
        function resetPrimal(obj)
            % erases existing data in primal variable storage
            obj.psiS = []; obj.thetaO = []; obj.thetaOB = [];
            obj.thetaS = []; obj.FS0 = []; obj.thetaBC1 = [];
            obj.thetaBC1 = []; obj.thetaBC3 = [];
            obj.thetaIC = []; obj.thetaT = []; obj.FT0 = [];
            obj.thetaC = []; obj.FC0 = []; obj.thetaChat = [];
            obj.theta = {}; obj.thetaStore = []; obj.ubar = [];
            obj.wbar = []; obj.ur = {}; obj.uz = {};
            obj.beta = []; obj.eta = []; obj.PsiCnm = []; obj.GCnm = [];
            obj.RD = []; obj.RFD = []; obj.RC3 = []; obj.RC4 = [];
            obj.eta = []; obj.beta = []; obj.RStatic = []; obj.cGet = [];
            obj.RC3Static = []; obj.RC4Static = []; obj.etaStatic = []; 
            obj.rhoH = []; obj.FH = [];
            obj.betaH = []; obj.etaH = []; obj.CnmH = [];
            obj.thetaH = []; obj.FHhatCn = []; obj.betaFH = [];   
            obj.thetaCh = []; obj.g1 = {}; obj.g3 = {};
            obj.g2 = {}; obj.g4 = {}; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % temperature and velocity simulation functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function simulateStorageTheta(obj, profiling, ...
                                         verify_conservation, hold_start)
            % simulates temperature in storage bin for arbitrary sequence
            % of charging/holding/discharging which is defined in the
            % FoMode property
            if nargin < 2, profiling = 0; end
            if nargin < 3, verify_conservation = 0; end
            if nargin < 4, hold_start = []; end
            % initialize temperature matrices
            ns = length(obj.zhat); nl = length(obj.zbar);
            ms = length(obj.rhat); ml = length(obj.rbar);
            mt = length(obj.rtop); nc = length(obj.zcenter);
            if isempty(obj.thetaS), obj.thetaS = ones(nl, ml); end
            if isempty(obj.thetaT), obj.thetaT = ones(ns, mt); end
            if isempty(obj.thetaC)
                obj.thetaC = ones(nc, ms);
                computeThetaChat(obj);
            end
            obj.FS0 = obj.thetaS;
            nch = length(obj.zbarH); mch = length(obj.rbarH);
            IC = ones(nch, mch); zIC = obj.zbarH; rIC = obj.rbarH;
            obj.Tp = obj.theta2TK(mean(IC(end, :)));
            % initialize wall systems
            initializeBaseSys(obj);
%             initializeTopSys(obj);
%             matchWallBoundary(obj, zIC, IC(:, end));
%             computeThetaW(obj, obj.Fo(1));
            initializeWallSys(obj);
%             initializeRoofSys(obj);
            initializeER(obj);
            % iterate through each time step
            kStart = 2;
            if profiling, profile on; end
            if ~isempty(hold_start)
                obj.theta = hold_start;
                IC = full(obj.theta{end, 1});
                zIC = obj.theta{end, 2}; rIC = obj.theta{end, 3};
                kStart = 2;
            end                   
            for kf = kStart:length(obj.Fo)
                obj.currentStatus = printStatus(obj, IC);
                obj.FoNow = obj.Fo(kf);
                obj.FoModeNow = obj.FoMode{1, kf};
                % normalize time step
                kn_ = mod(kf-1, obj.ls) + 1;
                switch obj.FoModeNow
                    case 'H'
                        % set initial condition
                        if kf == 1, prev = 1; else, prev = kf - 1; end
                        if obj.FoMode{1, prev} ~= 'H'
                            matchHoldIC(obj, IC, zIC, rIC, ...
                                                  verify_conservation);
                        end     
                        obj.Bi5 = obj.Bi5H;
                        % compute stagnant holding solution for entire
                        % holding sequence
                        computeThetaH(obj, obj.FoNow);
                        % store temperature values in solution matrix
                        patchThetaH(obj, kf, verify_conservation);
                        % update wall and base boundaries
                        ub = mean(mean(obj.theta{end, 1}(1, :)))*ones(10, 1);
                        uw = mean(mean(obj.theta{end, 1}(:, end)))*ones(10, 1);
                        obj.Tp = obj.theta2TK(mean(IC(end, :)));
%                         obj.Tp = obj.theta2TK(1);
                        Fo_ = linspace(0, obj.df, 10);
                        computeBaseSys(obj, ub, Fo_);                       
                        computeWallSys(obj, uw, Fo_); 
                        computeER(obj, Fo_);
                    case 'C'
                        % set initial condition
                        if kf == 1, prev = 1; else, prev = kf - 1; end
                        if obj.FoMode{1, prev} ~= 'C'
                            matchChargeIC(obj, IC, zIC, rIC, ...
                                                  verify_conservation);
                        end            
                        obj.Bi5 = obj.Bi5C;
                        % compute charge temperature values for current
                        % time step
                        computeThetaH(obj, obj.FoNow);
                        % store temperature values in solution matrix
                        patchThetaC(obj, kf, verify_conservation);
                        % update domain according to mass accounting
                        if kf ~= length(obj.Fo)
                            if obj.FoMode{1, kf+1} == 'C'                                              
                                updateDomainC(obj);
                            end  
                        end
                        % update wall and base boundaries
                        ub = mean(mean(obj.theta{end, 1}(1, :)))*ones(10, 1);
                        uw = mean(mean(obj.theta{end, 1}(:, end)))*ones(10, 1);
                        obj.Tp = obj.theta2TK(mean(IC(end, :)));
%                         obj.Tp = obj.theta2TK(1);
                        Fo_ = linspace(0, obj.df, 10);
                        computeBaseSys(obj, ub, Fo_);                       
                        computeWallSys(obj, uw, Fo_);
                        computeER(obj, Fo_);
                    case 'D'
                        % set initial condition
                        if kf == 1, prev = 1; else, prev = kf - 1; end
                        if obj.FoMode{1, prev} ~= 'D'
                            matchDischargeIC(obj, IC, zIC, rIC, ...
                                                  verify_conservation);
                        end 
                        obj.Bi5 = obj.Bi5D;                        
                        % compute temperature in top boundary
                        iterateThetaT(obj);
                        recordTopLoss(obj);
                        % compute temperature in center boundary
                        iterateThetaC(obj);
                        % compute the stagnant boundary offset temperature
                        computeThetaI(obj);
                        % compute temperature in stagnant region
                        iterateThetaS(obj);
                        % store temperature values in solution matrix
                        patchTheta(obj, kf, verify_conservation);
                        % update domain according to mass accounting
                        if kf ~= length(obj.Fo)
                            if obj.FoMode{1, kf+1} == 'D'
                                updateDomain(obj);
                            end
                        end
                        % update wall and base boundaries
                        ub = mean(mean(obj.theta{end, 1}(1, :)))*ones(10, 1);
                        uw = mean(mean(obj.theta{end, 1}(:, end)))*ones(10, 1);
                        obj.Tp = obj.theta2TK(mean(IC(end, :)));
%                         obj.Tp = obj.theta2TK(1);
                        Fo_ = linspace(0, obj.df, 10);
                        computeBaseSys(obj, ub, Fo_);                       
                        computeWallSys(obj, uw, Fo_);
                        computeER(obj, Fo_);
                end
                IC = full(obj.theta{kn_, 1});
                zIC = obj.theta{kn_, 2}; rIC = obj.theta{kn_, 3};
                reInitObj(obj);
                % compute remaining domain temperatures
%                 ut = obj.thetaA*ones(2, 1);
%                 computeTopSys(obj, ut, Fo_);
%                 computeThetaA(obj, IC, zIC, rIC);
%                 matchWallBoundary(obj, zIC, IC(:, end));
%                 computeThetaW(obj, obj.Fo(kf));                
                % hold if close to empty
                if obj.ztop < 0.1 && kf ~= length(obj.Fo)
                    if obj.FoMode{1, kf+1} == 'D'
                        obj.FoMode{1, kf+1} = 'H';
                    end                    
                end
                if obj.ztop > 0.95 && kf ~= length(obj.Fo)
                    if obj.FoMode{1, kf+1} == 'C'
                        obj.FoMode{1, kf+1} = 'H';
                    end                    
                end
            end  
            if profiling, profsave; end
        end
        function simulateCliffsModel(obj, zTop)
            % simulates the isolated wall model for a given time sequence,
            % obj.Fo, initial temperature distribution, obj.rhoW, and transient
            % boundary condition, TBC (obj.gW1)            
            obj.Bi5 = obj.Bi5H; 
            obj.zbarH = obj.zbarW;
            for k_ = 1:length(obj.Fo)
                obj.FoNow = obj.Fo(k_);
                obj.ztop = zTop(k_);
                if zTop(k_) == 0
                        Fo_ = obj.Fo(k_);
                        % compute wall solution                   
                        obj.gW1 = obj.Biw1*obj.thetaA ...
                                            *ones(length(obj.zbarW), 1);
                        obj.gW1 = obj.Biw1*obj.thetaA;
                        computeThetaW_new(obj, Fo_);
                        thetaP = obj.thetaA*ones(length(obj.zbarH), ...
                                                        length(obj.rbarH));
                else
                        % determine sequential holding sequence
                        Fo_ = obj.Fo(k_);                      
                        % compute wall solution                   
                        obj.gW1 = obj.Biw1*ones(length(obj.zbarW), 1);
%                         obj.gW1 = 1;
                        computeThetaW_new(obj, Fo_);  
                        thetaP = ones(length(obj.zbarH), length(obj.rbarH));
                end
                patchThetaK(obj, k_, thetaP);
            end
        end
        function simulateKevinsModel(obj, zTop)
            % simulates the isolated wall model for a given time sequence,
            % obj.Fo, initial temperature distribution, obj.rhoW, and transient
            % boundary condition, TBC (obj.gW1)            
            obj.Bi5 = obj.Bi5H; 
            obj.zbarH = obj.zbarW;
            for k_ = 1:length(obj.Fo)
                obj.FoNow = obj.Fo(k_);
%                 if k_ == 1
%                     obj.gW1 = 0.1*ones(length(obj.zbarW), 1);
%                     buildWall(obj); computeFW(obj); 
%                     obj.thetaW = obj.FW; obj.rhoW = obj.thetaW;
%                 end
                if zTop(k_) == 0                    
                    obj.rhoH = ones(length(obj.zbarH), length(obj.rbarH));
                    % compute wall solution                   
                    obj.gW1 = obj.thetaA*ones(length(obj.zbarW), 1);
                    computeThetaW(obj);
                    thetaP = obj.thetaA*ones(length(obj.zbarH), ...
                                                       length(obj.rbarH));
                else
                    % compute stagnant holding solution
%                     computeThetaH(obj, obj.FoNow);
                    % compute wall solution                   
%                     obj.gW1 = obj.thetaH(:, end);
                    obj.gW1 = ones(length(obj.zbarW), 1);
                    computeThetaW(obj, obj.df);  
%                     thetaP = obj.thetaH;
                    thetaP = ones(length(obj.zbarH), length(obj.rbarH));
                end
                patchThetaK(obj, k_, thetaP);
            end
        end        
        function matchWallBoundary(obj, z_, theta_)
            % generates the inner wall prescribed temperature BC with the
            % current air and particle temperatures
            [~, h_] = min(abs(z_(end) - obj.zbarW));
            gWs = interp1(z_, theta_, obj.zbarW(1:h_), 'makima');
%             gWs = mean(gWs);
            obj.gW1 = [gWs'.*ones(h_, 1); ...
                      obj.thetaA*ones(length(obj.zbarW) - h_, 1)];
            % add ramp function to smooth response
%             obj.gW1 = obj.gW1*(1 - exp(-obj.FoNow/obj.tauW1));
        end
        function matchHoldIC(obj, IC, zIC, rIC, verify_conservation)
            if nargin < 4, verify_conservation = 0; end
            obj.nzH = ceil(obj.nzH0*obj.ztop);
            [obj.zbarH, obj.dzH] = ...
                              nodeGen(obj, [0, obj.ztop], obj.nzH);
            [R, Z] = meshgrid(rIC, zIC);
            [Rq, Zq] = meshgrid(obj.rbarH, obj.zbarH);
            obj.rhoH = interp2(R, Z, IC, Rq, Zq, 'makima');
            % match dimensions for heat loss on top surface
%             obj.qTopP{1} = interp1(rIC, obj.qTopP{1}, obj.rbarH, 'makima');
%             obj.qTopP{2} = interp1(rIC, obj.qTopP{2}, obj.rbarH, 'makima');
            if verify_conservation
                % set matching temperature solution for energy equation
                obj.thetaMatchH{1} = obj.thetaMatchD{2}; 
            end
        end
        function matchChargeIC(obj, IC, zIC, rIC, verify_conservation)
            if nargin < 4, verify_conservation = 0; end
            obj.nzH = ceil(obj.nzH0*obj.ztop);
            [obj.zbarH, obj.dzH] = ...
                              nodeGen(obj, [0, obj.ztop], obj.nzH);
            [R, Z] = meshgrid(rIC, zIC);
            [Rq, Zq] = meshgrid(obj.rbarH, obj.zbarH);
            obj.rhoH = interp2(R, Z, IC, Rq, Zq, 'makima');
            obj.rhoH = ones(size(obj.rhoH));
%             % match dimensions for heat loss on top surface
%             obj.qTopP{1} = interp1(rIC, obj.qTopP{1}, obj.rbarH, 'makima');
%             obj.qTopP{2} = interp1(rIC, obj.qTopP{2}, obj.rbarH, 'makima');
            if verify_conservation
                % set matching temperature solution for energy equation
                obj.thetaMatchH{1} = obj.thetaMatchD{2}; 
            end
        end
        function matchDischargeIC(obj, IC, zIC, rIC, verify_conservation)
            if nargin < 4, verify_conservation = 0; end
            obj.ztop = obj.ztop - obj.h;
            obj.HNow = obj.ztop*obj.H;
            obj.nzbar = ceil(obj.nzbar0*obj.ztop);
            obj.nzc = ceil(obj.nzc0*obj.ztop);
            [obj.zbar, obj.dzbar] = ...
                            nodeGen(obj, [0, obj.ztop], obj.nzbar);
            [obj.zcenter, obj.dzc] = ...
                            nodeGen(obj, [0, obj.ztop], obj.nzc);
            [~, i] = min(abs(obj.zbar(end) - zIC));
            [~, j] = min(abs(obj.rbar(1) - rIC));
            % match dimensions for stagnant region            
            [R, Z] = meshgrid(rIC(j:end), zIC(1:i));
            [Rq, Zq] = meshgrid(obj.rbar, obj.zbar);
            obj.thetaS = interp2(R, Z, IC(1:i, j:end), Rq, Zq, 'makima');
            obj.FS0 = obj.thetaS;
            % match dimensions for center channel           
            [R, Z] = meshgrid(rIC(1:j), zIC(1:i));
            [Rq, Zq] = meshgrid(obj.rhat, obj.zcenter);
            obj.thetaC = interp2(R, Z, IC(1:i, 1:j), Rq, Zq, 'makima');
            % match dimensions for top flow surface
            [R, Z] = meshgrid(rIC(j:end), zIC(i:end) - zIC(i));
            [Rq, Zq] = meshgrid(obj.rtop, obj.zhat);
            obj.thetaT = interp2(R, Z, IC(i:end, j:end), Rq, Zq, 'makima');
            yt = obj.T2S('top'); yc = obj.C2S('top');          
            obj.g1 = obj.Bi1*(ones(size(yt)) - obj.thetaI); 
            obj.g3 = obj.Bi3*(ones(size(yc)) - obj.thetaI); 
            % match dimensions for heat loss on top surface
%             obj.qTopP{1} = interp1(rIC, obj.qTopP{1}, ...
%                                  [obj.rhat(1:end-1), obj.rbar], 'makima');
%             obj.qTopP{2} = interp1(rIC, obj.qTopP{2}, ...
%                                  [obj.rhat(1:end-1), obj.rbar], 'makima');
            if verify_conservation
                % set matching temperature solution for energy equation
                obj.thetaMatchD{1} = obj.thetaMatchH{2}; 
            end
        end
        function simulateChargeTheta(obj, profiling)
            % main function for computing the 2D temperature distribution 
            % in charging mode
            if nargin < 2, profiling = 0; end
            % set initial charging condition
            obj.ztop = obj.z0C;
            reInitObj(obj);
            % initialize temperature matrices
            n = length(obj.zbarH); m = length(obj.rbarH);
            if isempty(obj.rhoH), obj.rhoH = ones(n, m); end
            obj.thetaCh = obj.rhoH;
            % iterate through each time step
            if profiling, profile on; end
            for kf = 1:length(obj.Fo)
                obj.FoNow = obj.Fo(kf); 
                % patch and save in theta cell array
                patchThetaC(obj, kf);
                % compute distribution for current domain
                iterateThetaCh(obj);
                % update domain according to mass accounting
                updateDomainC(obj);                
            end  
            if profiling, profsave; end            
        end
        function simulateHoldTheta(obj, profiling)
            % main function for computing the transient 2D temperature 
            % profile throughout a storage holding period  
            if nargin < 2, profiling = 0; end
            % initialize temperature matrices
            n = length(obj.zbarH); m = length(obj.rbarH);
            if isempty(obj.rhoH), obj.rhoH = ones(n, m); end            
            if profiling, profile on; end
            % segment time array based on storage allocation
            nk = ceil(length(obj.Fo)/obj.ls);
            for k_ = 1:nk
                try % set full time segment
                    ksplit = (k_-1)*obj.ls+1:k_*obj.ls;
                    Fo_ = obj.Fo(ksplit);
                catch % set remaining time segment
                    ksplit = (k_-1)*obj.ls+1:length(obj.Fo);
                    Fo_ = obj.Fo(ksplit);
                end
                % compute analytical solution for current time segment
                computeThetaH(obj, Fo_);
                patchThetaH(obj, ksplit);
                saveTheta(obj, ksplit(end));                
            end
            if profiling, profsave; end
        end
        function simulateDischargeTheta(obj, storage, profiling, ...
                                         verify_conservation, sensitivity)
            % main function for coupling stagnant temperature solution with
            % the boundary solutions to simulate the temperature
            % propagation in the full domain for storage discharge 
            if nargin < 2, storage = 'full'; end
            if nargin < 3, profiling = 0; end
            if nargin < 4, verify_conservation = 0; end
            if nargin < 5, sensitivity = 0; end
            % initialize temperature matrices
            ns = length(obj.zhat); nl = length(obj.zbar);
            ms = length(obj.rhat); ml = length(obj.rbar);
            mt = length(obj.rtop); nc = length(obj.zcenter);
            if isempty(obj.thetaS), obj.thetaS = ones(nl, ml); end
            if isempty(obj.thetaT), obj.thetaT = ones(ns, mt); end
            if isempty(obj.thetaC)
                obj.thetaC = ones(nc, ms);
                computeThetaChat(obj);
            end
            obj.FS0 = obj.thetaS;
            % iterate through each time step
            if profiling, profile on; end
            for kf = 1:length(obj.Fo)
                obj.FoNow = obj.Fo(kf); 
                % store temperature values in solution matrix
                switch storage
                    case 'full'
                        patchTheta(obj, kf, verify_conservation,...
                                                          sensitivity);    
                    case 'point'
                        savePoints(obj);
                end
                % compute temperature in top boundary
                iterateThetaT(obj);
                recordTopLoss(obj);
                % compute temperature in center boundary
                iterateThetaC(obj, sensitivity);   
                % compute temperature in stagnant regiontheta
                iterateThetaS(obj);   
                % update domain according to mass accounting
                updateDomain(obj);                
            end  
            if profiling, profsave; end
            obj.ztop = 1;
            reInitObj(obj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conduction temperature solution for holding and charging mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeThetaH(obj, Fo_)
            % compute filter function solution
            computeFH(obj); % needs to be before eigenvalue comp.
            % compute eigenvalues and coefficients if not yet populated
            computeBetaH(obj);
            computeEtaH(obj);
            computeCnmH(obj);
            if nargin < 2, Fo_ = obj.Fo; end
            n = length(obj.zbarH); m = length(obj.rbarH); k_ = length(Fo_);             
            % compute Fourier series for filtered solution
            obj.thetaH = sparse(zeros(n*k_, m*k_));
            for i = 1:length(obj.CnmH)
                Z = XHn(obj, obj.zbarH, obj.betaH(i));
                R = XHm(obj, obj.rbarH, obj.etaH(i));
                F = Xt(obj, Fo_, obj.betaH(i), obj.etaH(i));
                obj.thetaH = obj.thetaH + ...
                             sparse(kron(obj.CnmH(i)*diag(F), Z'*R));               
            end             
            % shift back with filtering function to get original solution
            obj.thetaH = obj.thetaH - ...
                                sparse(kron(eye(length(Fo_)), obj.FH));
            obj.rhoH = obj.thetaH((k_-1)*n+1:end, (k_-1)*m+1:end);
        end
        function computeFH(obj)
            % computes the filtering function for the bin holding solution           
            if isempty(obj.g2), initializeBaseSys(obj); end
            if isempty(obj.g4), initializeWallSys(obj); end
            computeFHhatCn(obj);    
            % compute FHhat
            n = length(obj.zbarH); m = length(obj.rbarH);
            % first s.s. BVP filter
            FH1 = zeros(n, m);
            for i = 1:length(obj.FHhatCn)
                Z = XHn(obj, obj.zbarH, obj.betaFH(i));
                R = XFHm(obj, obj.rbarH, obj.betaFH(i));
                FH1 = FH1 + obj.FHhatCn(i)*Z'*R;               
            end
            % second s.s. BVP filter
            FH2 = repmat(obj.g2*(obj.zbarH' - obj.ztop - 1/obj.Bi5), 1, m);
            % full filter superposition
            obj.FH = FH1 + FH2 - obj.thetaA;
        end
        function xi = XHn(~, zbar_, beta_)
            % eigenfunction for z-dimension BVP
            xi = cos(beta_*zbar_);
        end
        function ni = NHn(obj, beta_)
            % norm for z-dimension BVP
            ni = 0.5*(obj.ztop*(beta_.^2 + obj.Bi5^2) + obj.Bi5)./ ...
                      (beta_.^2 + obj.Bi5^2);
        end
        function zi = ZHn(obj, beta_)
            % transendental function for beta values
            zi = beta_.*sin(beta_*obj.ztop) - obj.Bi5*cos(beta_*obj.ztop);
        end 
        function xi = XHm(~, rbar_, eta_)
            % eigenfunction for r-dimension BVP
            xi = besselj(0, eta_*rbar_);
        end
        function ni = NHm(obj, eta_)
            % norm for r-dimension BVP
            if eta_ == 0
                ni = 0.5*obj.b^2;
            else
                ni = 0.5*obj.b^2*besselj(0, eta_*obj.b).^2;
            end
        end
        function zi = ZHm(obj, eta_)
            % transendental function for eta values
            zi = besselj(1, eta_*obj.b);
        end
        function xi = XFHn(obj, zbar_, eta_)
            % z-BVP solution for holding filter 2
            xi = eta_*(sinh(eta_*obj.ztop)*sinh(eta_*zbar_) - ...
                       cosh(eta_*obj.ztop)*cosh(eta_*zbar_)) + ...
                 obj.Bi5*(sinh(eta_*obj.ztop)*sinh(eta_*zbar_) - ...
                         (eta_*obj.ztop)*cosh(eta_*zbar_));
        end
        function xi = XFHm(~, rbar_, beta_)
            % r-BVP solution for holding filter 1
            xi = besseli(0, beta_*rbar_);
        end
        function c = FourierCoefficientThetaHhat(obj, beta_, eta_)
            % Fourier coefficients
            C_num = computeThetaHhatIntegral(obj, obj.zbarH, obj.rbarH, ...
                                                            beta_, eta_);
            c = C_num/(NHm(obj, eta_)*NHn(obj, beta_));
        end
        function c = FourierCoefficientFH1(obj, beta_)
            % Fourier coefficients
            C_num = mean(mean(obj.g4))*sin(beta_*obj.ztop)./ ...
                                      (beta_.^2.*besseli(1, beta_*obj.b));
            c = C_num/NHn(obj, beta_);  
        end
        function x = computeThetaHhatIntegral(obj, z_, r_, beta_, eta_)
            % computes integral from numerator of Fourier coefficient for
            % the filtered thetaH solution
            if isempty(obj.rhoH)
                obj.rhoH = ones(length(z_), length(r_));
            end
            if isempty(obj.FH), computeFH(obj); end
            fr = simpsonIntegrator(obj, r_);
            fz = simpsonIntegrator(obj, z_);
            Ir = (XHm(obj, r_, eta_).*r_.*(obj.rhoH + obj.FH))*fr';
            x = (XHn(obj, z_, beta_).*Ir')*fz';
        end                
        % computation of beta values
        function computeBetaH(obj, ShowPlot)
            interval = linspace(0, obj.pfH, obj.pH); % interval/spacing 
                                                     % of root calculation
            rn = NaN*ones(obj.pH, 1);                % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.pH
                rn(i) = fzero(@(beta_)ZHn(obj, beta_), interval(i), options);
            end
            obj.betaH = rn(diff(rn)>1e-10);        % only keep unique roots
            obj.betaH = obj.betaH(obj.betaH > 1e-10); % remove zeros  
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(interval, ZHn(obj, interval), '-k');
                hold on
                plot(obj.betaH, ZHn(obj, obj.betaH), '.r');
                legend('$ZH_n(\hat{\beta})$', '$\beta_n$', 'interpreter', ...
                    'latex', 'FontSize', 14, 'NumColumns', 2);
                xlabel('$\hat{\beta}$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylabel('$ZH_n(\hat{\beta})$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                xlim([0, obj.pfH]);
                hold off 
            end                   
        end 
        function computeBetaFH(obj, ShowPlot)
            interval = linspace(0, obj.pfFH, obj.pFH);   % interval/spacing 
                                                     % of root calculation
            rm = NaN*ones(obj.pFH, 1);              % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.pFH
                rm(i) = fzero(@(beta_) ZHn(obj, beta_), interval(i), options);
            end
            obj.betaFH = rm(diff(rm)>1e-10);     % only keep unique roots
            obj.betaFH = obj.betaFH(obj.betaFH > 1e-10);   % remove zeros
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', ...
                    'Position', [0 0 0.4 0.25], 'Visible', 'on');
                plot(interval, ZHn(obj, interval), '-k');
                hold on
                plot(obj.betaFH, ZHn(obj, obj.betaFH), '.r');
                legend('$ZH_n(\hat{\eta})$', '$\beta_n$', 'interpreter',...
                    'latex', 'FontSize', 14, 'NumColumns', 2);
                xlabel('$\hat{\beta}$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$ZH_n(\hat{\beta})$', 'interpreter', 'latex',...
                    'FontSize', 14);
                xlim([0, obj.pfFH]);
                hold off             
            end    
        end 
        % computation of eta values
        function computeEtaH(obj, ShowPlot)
            interval = linspace(0, obj.qfH, obj.qH);   % interval/spacing 
                                                     % of root calculation
            rm = NaN*ones(obj.qH, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.qH
                    rm(i) = fzero(@(eta_) ZHm(obj, eta_), ...
                                                    interval(i), options);
            end
            obj.etaH = rm(diff(rm) > 1e-10);       % only keep unique roots
            obj.etaH(1) = 0;
            obj.etaHStatic = obj.etaH;
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(interval, ZHm(obj, interval), '-k');
                hold on
                plot(obj.etaH, ZHm(obj, obj.etaH), '.r');
                legend('$ZH_m(\hat{\eta})$', '$\eta_m$', 'interpreter',...
                    'latex', 'FontSize', 14, 'NumColumns', 2);
                xlabel('$\hat{\eta}$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$ZH_m(\hat{\eta})$', 'interpreter', 'latex',...
                    'FontSize', 14);
                xlim([0, obj.qfH]);
                hold off 
            end    
        end    
        % computation of fourier coefficients
        function computeFHhatCn(obj, ShowPlot)
            % compute beta and eta values if not already populated
            if nargin > 1 && ShowPlot
                computeBetaFH(obj, true);
            else
                computeBetaFH(obj);
            end
            if length(obj.betaFH) < obj.niFH
                obj.niFH = length(obj.betaFH);
            end
            CnTemp = NaN*ones(1, obj.niFH);
            betaTemp = NaN*ones(1, obj.niFH);            
            for j = 1:obj.niFH
                CnTemp(j) = FourierCoefficientFH1(obj, obj.betaFH(j));
                betaTemp(j) = obj.betaFH(j);
            end
            cGet_ = abs(CnTemp) > obj.climFH;
            obj.FHhatCn = CnTemp(cGet_);
            obj.betaFH = betaTemp(cGet_);
            cGetIdx = find(cGet_);                      
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(cGetIdx, obj.FHhatCn, 'rx');
                hold on
                plot(1:length(CnTemp(:)), CnTemp(:), '-k');
                legend('Used', 'Computed', 'interpreter', 'latex', ...
                    'FontSize', 14, 'NumColumns', 2);
                xlabel('$n$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$C_{n}(\beta_n)$', 'interpreter', ...
                    'latex', 'FontSize', 14);
                xlim([0, obj.niFH]);
                hold off 
            end               
        end
        function computeCnmH(obj, ShowPlot)
            % compute beta and eta values if not already populated
            if isempty(obj.betaH)
                if nargin > 1 && ShowPlot
                    computeBetaH(obj, true);
                else
                    computeBetaH(obj);
                end
            end
            if isempty(obj.etaHStatic)
                if nargin > 1 && ShowPlot
                    computeEtaH(obj, true);
                else
                    computeEtaH(obj);
                end
            end
            if length(obj.betaH) < obj.niH 
                ni_ = length(obj.betaH);
            else
                ni_ = obj.niH;
            end
            CnmTemp = NaN*ones(ni_, obj.miH);
            betaTemp = NaN*ones(ni_, obj.miH);
            etaTemp = NaN*ones(ni_, obj.miH);            
            for i = 1:ni_
                for j = 1:obj.miH
                    CnmTemp(i, j) = FourierCoefficientThetaHhat(obj, ...
                        obj.betaH(i), obj.etaHStatic(j));
                    betaTemp(i, j) = obj.betaH(i);
                    etaTemp(i, j) = obj.etaHStatic(j);
                end
            end
            cGet_ = abs(CnmTemp) > obj.climH;
            obj.CnmH = CnmTemp(cGet_);
            obj.betaH = betaTemp(cGet_);
            obj.etaH = etaTemp(cGet_);
            cGetIdx = find(cGet_);  
%             w = obj.betaH*pi/max(obj.betaH); 
%             sigmaN = sin(w)./w;
            sigmaN = ones(length(obj.betaH), 1);
            obj.CnmH = sigmaN.*obj.CnmH;
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(cGetIdx, obj.CnmH, 'rx');
                hold on
                plot(1:length(CnmTemp(:)), CnmTemp(:), '-k');
                legend('Used', 'Computed', 'interpreter', 'latex', ...
                    'FontSize', 14, 'NumColumns', 2);
                xlabel('$nm$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$C_{nm}(\beta_n, \eta_m)$', 'interpreter', ...
                    'latex', 'FontSize', 14);
                xlim([0, obj.niH*obj.miH]);
                hold off 
            end               
        end         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conduction temperature solution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function iterateThetaS(obj)      
            % compute auxilary problem solution 
            % compute homogeneous filters
            computeFD(obj);
            % eigenvalues and Fourier coeffiecients
            computeGCnm(obj);
            % initial condition contribution
            computeThetaIC(obj);
            % compute boundary condition contributions
            computeThetaBCNow(obj);            
            % sum individual solution terms
            obj.thetaS = obj.thetaIC + ...
                                   obj.thetaBC1 + obj.thetaBC3 - obj.FD;                                       
            obj.FS0 = obj.thetaS;
        end          
        % z-dimension temperature BVP
        function xi = Xn(~, zbar_, beta_)
            % eigenfunction for z-dimension BVP
            xi = cos(beta_*zbar_);
        end
        function ni = Nn(obj, beta_)
            % norm for z-dimension BVP
            ni = 0.5*(obj.ztop*(beta_.^2 + obj.Bi1^2) + obj.Bi1)./ ...
                                                (beta_.^2 + obj.Bi1^2);
        end
        function zi = Zn(obj, beta_)
            % transendental function for beta values
            zi = beta_.*sin(obj.ztop*beta_) - obj.Bi1*cos(obj.ztop*beta_);
        end     
        function xi = XFn(obj, zbar_, eta_)
            % z-BVP eigenfunction for filter 2
            xi = eta_.*(sinh(eta_*obj.ztop).*sinh(eta_*zbar_) - ...
                        cosh(eta_*obj.ztop).*cosh(eta_*zbar_)) + ...
               obj.Bi1*(cosh(eta_*obj.ztop).*sinh(eta_*zbar_) - ...
                        sinh(eta_*obj.ztop).*cosh(eta_*zbar_));
        end
        % r-dimension temperature BVP
        function si = S(obj, eta_)
            si = -eta_.*bessely(1, eta_*obj.b) + ...
                obj.Bi4*bessely(0, eta_*obj.b);
        end
        function ui = U(obj, eta_)
            ui = -eta_.*besselj(1, eta_*obj.a0) - ...
                obj.Bi3*besselj(0, eta_*obj.a0);
        end
        function vi = V(obj, eta_)
            vi = -eta_.*besselj(1, eta_*obj.b) + ...
                obj.Bi4*besselj(0, eta_*obj.b);
        end
        function wi = W(obj, eta_)
            wi = -eta_.*bessely(1, eta_*obj.a0) - ...
                obj.Bi3*bessely(0, eta_*obj.a0);
        end
        function b1i = B1(obj, eta_)
            b1i = obj.Bi3^2 + eta_.^2;
        end
        function b2i = B2(obj, eta_) 
            b2i = obj.Bi4^2 + eta_.^2;
        end
        function xi = Xm(obj, rbar_, eta_)
            % eigenfunction for r-dimension BVP
            xi = besselj(1, eta_*obj.b)*bessely(0, eta_*rbar_) - ...
                 besselj(0, eta_*rbar_)*bessely(1, eta_*obj.b);
        end
        function ni = Nm(obj, eta_)
            % norm for r-dimension BVP
            ni = (2/pi^2)*(U(obj, eta_).^2 - B1(obj, eta_).* ...
                   besselj(1, eta_*obj.b).^2)./(eta_.^2.*U(obj, eta_).^2);
        end
        function zi = Zm(obj, eta_)
            % transendental function for eta values
            zi = W(obj, eta_).*besselj(1, eta_*obj.b) - ...
                 U(obj, eta_).*bessely(1, eta_*obj.b);
        end   
        function [xi, xip] = XFm(obj, rbar_, beta_)
            % eigenfunction for r-component of filter 1
            xi = obj.Bi3*(besselk(0, beta_*rbar_).* ...
                          besseli(0, beta_*obj.a0) - ...
                          besselk(0, beta_*obj.a0).* ...
                          besseli(0, beta_*rbar_)) + ...
                 beta_*(besselk(0, beta_*rbar_).* ...
                        besseli(1, beta_*obj.a0) - ...
                        besselk(1, beta_*obj.a0).* ...
                        besseli(0, beta_*rbar_));
            xip = -beta_*(obj.Bi3*(besselk(1, beta_*rbar_).* ...
                                   besseli(0, beta_*obj.a0) + ...
                                   besselk(0, beta_*obj.a0).* ...
                                   besseli(1, beta_*rbar_)) + ...
                          beta_*(besselk(1, beta_*rbar_).* ...
                                 besseli(1, beta_*obj.a0) + ...
                                 besselk(1, beta_*obj.a0).* ...
                                 besseli(1, beta_*rbar_)));            
        end
        % combined homogeneous temperature solution           
        function c = FourierCoefficientG(obj, beta_, eta_, R_)
            % Fourier coefficients
            C_num = computeICIntegral(obj, obj.zbar, obj.rbar, beta_, R_);
            c = C_num/(Nn(obj, beta_)*Nm(obj, eta_));
        end   
        function c = FourierCoefficientFD1(obj, beta_)
            % Fourier coefficients
            [~, Rpb] = XFm(obj, obj.b, beta_);
            C_num = mean(mean(obj.g4))*sin(beta_*obj.ztop)/(beta_*Rpb);
            c = C_num/Nn(obj, beta_);  
        end
        function c = FourierCoefficientFD2(obj, eta_)
            % Fourier coefficients
            C_num = obj.g2*(besselj(1, eta_*obj.b).* ...
                                   (obj.b*bessely(1, eta_*obj.b) - ...
                                    obj.a0*bessely(1, eta_*obj.a0)) - ...
                            bessely(1, eta_*obj.b).* ...
                                   (obj.b*besselj(1, eta_*obj.b) - ...
                                    obj.a0*besselj(1, eta_*obj.a0)))/ ...
                    (eta_.^2.*(eta_.*sinh(eta_*obj.ztop) + ...
                               obj.Bi1*cosh(eta_*obj.ztop)));
            c = C_num/Nm(obj, eta_);  
        end
        function computeThetaIC(obj)
            % compute eigenvalues and coefficients if not yet populated         
            n = length(obj.zbar); m = length(obj.rbar);
            % compute for current time step
            obj.thetaIC = zeros(n, m);
            obj.thetaICDot = zeros(n, m);
            for i = 1:length(obj.GCnm)
                C = obj.GCnm(i);
                F = Xt(obj, obj.df, obj.beta(i), obj.eta(i));
                Z = Xn(obj, obj.zbar, obj.beta(i));          
                obj.thetaIC = obj.thetaIC + C*F*Z'*obj.RD{i};
            end 
            obj.thetaIC = obj.thetaIC + obj.thetaI;
        end
        function x = computeICIntegral(obj, z_, r_, beta_, R_)
            % computes integral component of the initial condition
            % contribution in Green's formula
            if isempty(obj.FS0)
                obj.FS0 = ones(length(z_), length(r_));
            end
            if isempty(obj.FD), computeFD(obj); end
            fr = simpsonIntegrator(obj, r_);
            fz = simpsonIntegrator(obj, z_);
            Ir = (R_.*r_.*(obj.FS0 + obj.FD - obj.thetaI))*fr';
            x = (Xn(obj, z_, beta_).*Ir')*fz';
        end          
        function computeFD(obj)
            % computes the filtering function for the bin holding solution
            if isempty(obj.g2), initializeBaseSys(obj); end
            if isempty(obj.g4), initializeWallSys(obj); end
            computeBetaFD(obj);
            if isempty(obj.etaFDStatic), computeEtaFD(obj); end
            computeFDCn(obj);  
            computeFDCm(obj);
            % compute FHhat
            n = length(obj.zbar); m = length(obj.rbar);
            FD1 = zeros(n, m); FD2 = zeros(n, m);
            for i = 1:length(obj.FDCn)
                Z = Xn(obj, obj.zbar, obj.betaFD(i));
                R = XFm(obj, obj.rbar, obj.betaFD(i));
                FD1 = FD1 + obj.FDCn(i)*Z'*R;               
            end
            for i = 1:length(obj.FDCm)
                Z = XFn(obj, obj.zbar, obj.etaFD(i));
                FD2 = FD2 + obj.FDCm(i)*Z'*obj.RFD{i};               
            end
            obj.FD = FD1 + FD2;
        end
        function computeThetaBCNow(obj)
            % computes boundary integration value at current time step
            n = length(obj.zbar); m = length(obj.rbar);
            obj.thetaBC1 = zeros(n, m); obj.thetaBC3 = zeros(n, m);
            yt = obj.T2S('top'); yc = obj.C2S('top');
            fr = simpsonIntegrator(obj, obj.rbar);
            fz = simpsonIntegrator(obj, obj.zbar);
            % compute top and center boundary temperatures 
            g1_ = obj.Bi1*(yt - obj.thetaI);
            g3_ = obj.Bi3*(yc - obj.thetaI);
            if isempty(obj.g1)
                obj.g1 = obj.Bi1*(ones(size(g1_)) - obj.thetaI); 
            end
            if isempty(obj.g3)
                obj.g3 = obj.Bi3*(ones(size(g3_)) - obj.thetaI); 
            end
            for i = 1:length(obj.GCnm)
                Z = Xn(obj, obj.zbar, obj.beta(i))';
                C1 = Xn(obj, obj.ztop, obj.beta(i)) ...
                    /(Nn(obj, obj.beta(i))*Nm(obj, obj.eta(i)));
                C3 = obj.a0*obj.RC3(i) ...
                    /(Nn(obj, obj.beta(i))*Nm(obj, obj.eta(i)));
                IBC1_ = (obj.g1.*obj.RD{i}.*obj.rbar)*fr';
                IBC3_ = (obj.g3(1:n).*Xn(obj, obj.zbar, obj.beta(i)))*fz';
                IBC1 = (g1_.*obj.RD{i}.*obj.rbar)*fr';
                IBC3 = (g3_.*Xn(obj, obj.zbar, obj.beta(i)))*fz';
                p1 = (IBC1 - IBC1_)/obj.df; q1 = IBC1_;
                p3 = (IBC3 - IBC3_)/obj.df; q3 = IBC3_;
                lambda = sqrt(obj.beta(i)^2 + obj.eta(i)^2);
                F1 = q1*(1 - exp(-obj.df*lambda^2))/lambda^2 + ...
                     p1*(1 - (lambda^2*obj.df + 1)* ...
                      exp(-obj.df*lambda^2))/lambda^4;
                F3 = q3*(1 - exp(-obj.df*lambda^2))/lambda^2 + ...
                     p3*(1 - (lambda^2*obj.df + 1)* ...
                      exp(-obj.df*lambda^2))/lambda^4;
                ep1 = C1*F1*Z*obj.RD{i}; 
                ep3 = C3*F3*Z*obj.RD{i}; 
                obj.thetaBC1 = obj.thetaBC1 + ep1;
                obj.thetaBC3 = obj.thetaBC3 + ep3;
            end
            obj.g1 = g1_; obj.g3 = g3_;
        end
        % computation of beta values
        function computeBeta(obj, ShowPlot)
            interval = linspace(0, obj.pf, obj.p);   % interval/spacing 
                                                     % of root calculation
            rn = NaN*ones(obj.p, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.p
                rn(i) = fzero(@(beta_)Zn(obj, beta_), interval(i), options);
            end
            obj.beta = rn(diff(rn)>1e-10);         % only keep unique roots
            obj.beta = obj.beta(obj.beta > 1e-10); % remove zeros            
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(interval, Zn(obj, interval), '-k');
                hold on
                plot(obj.beta, Zn(obj, obj.beta), '.r');
                legend('$Z_n(\hat{\beta})$', '$\beta_n$', 'interpreter', ...
                    'latex', 'FontSize', 14, 'NumColumns', 2);
                xlabel('$\hat{\beta}$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylabel('$Z_n(\hat{\beta})$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                xlim([0, obj.pf]);
                hold off 
            end                   
        end       
        function computeBetaFD(obj, ShowPlot)
            interval = linspace(0, obj.pfFD, obj.pFD);   % interval/spacing 
                                                     % of root calculation
            rn = NaN*ones(obj.pFD, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.pFD
                rn(i) = fzero(@(beta_)Zn(obj, beta_), interval(i), options);
            end
            obj.betaFD = rn(diff(rn)>1e-10);         % only keep unique roots
            obj.betaFD = obj.betaFD(obj.betaFD > 1e-10); % remove zeros            
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(interval, Zn(obj, interval), '-k');
                hold on
                plot(obj.betaFD, Zn(obj, obj.betaFD), '.r');
                legend('$Z_n(\hat{\beta})$', '$\beta_n$', 'interpreter', ...
                    'latex', 'FontSize', 14, 'NumColumns', 2);
                xlabel('$\hat{\beta}$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylabel('$Z_n(\hat{\beta})$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                xlim([0, obj.pfFD]);
                hold off
            end                   
        end   
        % computation of eta values
        function computeEta(obj, ShowPlot)
            interval = linspace(1, obj.qf, obj.q);   % interval/spacing 
                                                     % of root calculation
            rm = NaN*ones(obj.q, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.q
                rm(i) = fzero(@(eta_) Zm(obj, eta_), interval(i), options);
            end
            obj.eta = rm(diff(rm)>1e-10);         % only keep unique roots
            obj.eta = obj.eta(obj.eta > 1e-10);   % remove zeros  
            % save eta and eta dependencies for future reference
            obj.etaStatic = obj.eta;
            obj.RStatic = cell(length(obj.eta), 1);
            obj.RC3Static = NaN*ones(length(obj.eta), 1);
            obj.RC4Static = NaN*ones(length(obj.eta), 1);
            for i = 1:length(obj.eta)            
                obj.RStatic{i} = Xm(obj, obj.rbar, obj.eta(i));
                obj.RC3Static(i) = Xm(obj, obj.a0, obj.eta(i));
                obj.RC4Static(i) = Xm(obj, obj.b, obj.eta(i));
            end
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(interval, Zm(obj, interval), '-k');
                hold on
                plot(obj.eta, Zm(obj, obj.eta), '.r');
                legend('$Z_m(\hat{\eta})$', '$\eta_m$', 'interpreter', ...
                    'latex', 'FontSize', 14, 'NumColumns', 2);
                xlabel('$\hat{\eta}$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$Z_m(\hat{\eta})$', 'interpreter', 'latex',...
                    'FontSize', 14);
                xlim([0, obj.qf]);
                hold off 
            end    
        end       
        function computeEtaFD(obj, ShowPlot)
            interval = linspace(1, obj.qfFD, obj.qFD);   % interval/spacing 
                                                     % of root calculation
            rm = NaN*ones(obj.qFD, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.qFD
                rm(i) = fzero(@(eta_) Zm(obj, eta_), interval(i), options);
            end
            obj.etaFD = rm(diff(rm)>1e-10);         % only keep unique roots
            obj.etaFD = obj.etaFD(obj.etaFD > 1e-10);   % remove zeros  
            % save eta and eta dependencies for future reference
            obj.etaFDStatic = obj.etaFD;
            obj.RFDStatic = cell(length(obj.etaFD), 1);
            for i = 1:length(obj.etaFD)            
                obj.RFDStatic{i} = Xm(obj, obj.rbar, obj.etaFD(i));
            end
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(interval, Zm(obj, interval), '-k');
                hold on
                plot(obj.etaFD, Zm(obj, obj.eta), '.r');
                legend('$Z_m(\hat{\eta})$', '$\eta_m$', 'interpreter',...
                    'latex', 'FontSize', 14, 'NumColumns', 2);
                xlabel('$\hat{\eta}$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$Z_m(\hat{\eta})$', 'interpreter', 'latex',...
                    'FontSize', 14);
                xlim([0, obj.qfFD]);
                hold off 
            end    
        end   
        % computation of fourier coefficients
        function computeGCnm(obj, ShowPlot)
            % compute new beta values
            if nargin > 1 && ShowPlot
                computeBeta(obj, true);
            else
                computeBeta(obj);
            end
            % compute eta values if not already populated
            if isempty(obj.eta)
                if nargin > 1 && ShowPlot
                    computeEta(obj, true);
                else
                    computeEta(obj);
                end     
            end
            if length(obj.beta) < obj.ni
                obj.ni = length(obj.beta);
                warning('beta eigenvalue set has been shortened')
            end
            if length(obj.eta) < obj.mi
                obj.mi = length(obj.eta);
                warning('eta eigenvalue set has been shortened')
            end
            CnmTemp = NaN*ones(obj.ni, obj.mi);
            betaTemp = NaN*ones(obj.ni, obj.mi);
            etaTemp = NaN*ones(obj.ni, obj.mi);            
            RC3Temp = NaN*ones(obj.ni, obj.mi);
            RC4Temp = NaN*ones(obj.ni, obj.mi);
            RTemp = cell(obj.ni, obj.mi);
            for i = 1:obj.ni
                for j = 1:obj.mi
                    CnmTemp(i, j) = FourierCoefficientG(obj, ...
                        obj.beta(i), obj.etaStatic(j), obj.RStatic{j});
                    betaTemp(i, j) = obj.beta(i);
                    etaTemp(i, j) = obj.etaStatic(j);
                    RTemp{i, j} = obj.RStatic{j};
                    RC3Temp(i, j) = obj.RC3Static(j);
                    RC4Temp(i, j) = obj.RC4Static(j);
                end
            end
            obj.cGet = abs(CnmTemp) > obj.clim;
            GCnm_ = CnmTemp(obj.cGet);
            obj.beta = betaTemp(obj.cGet);
%             w = obj.beta*pi/max(obj.beta); 
%             sigmaN = sin(w)./w;
            sigmaN = ones(length(obj.beta), 1);
            obj.GCnm = sigmaN.*GCnm_;
            obj.eta = etaTemp(obj.cGet);
            obj.RD = RTemp(obj.cGet);
            obj.RC3 = RC3Temp(obj.cGet);
            obj.RC4 = RC4Temp(obj.cGet);
            cGetIdx = find(obj.cGet);  
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', ...
                    'Position', [0 0 0.4 0.25], 'Visible', 'on');
                plot(cGetIdx, GCnm_, 'rx');
                hold on
                plot(1:length(CnmTemp(:)), CnmTemp(:), '-k');
                legend('Used', 'Computed', 'interpreter', 'latex', ...
                    'FontSize', 14, 'NumColumns', 2);
                xlabel('$nm$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$C_{nm}(\beta_n, \eta_m)$', 'interpreter', ...
                    'latex', 'FontSize', 14);
                xlim([0, obj.ni*obj.mi]);
                hold off 
            end               
        end
        function computeFDCn(obj, ShowPlot)
            % compute beta and eta values if not already populated
            if nargin > 1 && ShowPlot
                computeBetaFD(obj, true);
            else
                computeBetaFD(obj);
            end 
            if length(obj.betaFD) < obj.niFD
                obj.niFD = length(obj.betaFD);
                warning('beta eigenvalue set has been shortened')
            end
            CnTemp = NaN*ones(1, obj.niFD);
            betaTemp = NaN*ones(1, obj.niFD);            
            for j = 1:obj.niFD
                CnTemp(j) = FourierCoefficientFD1(obj, obj.betaFD(j));
                betaTemp(j) = obj.betaFD(j);
            end
            cGet_ = abs(CnTemp) > obj.climFD;
            obj.FDCn = CnTemp(cGet_);
            obj.betaFD = betaTemp(cGet_);
            cGetIdx = find(cGet_);                      
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(cGetIdx, obj.FDCn, 'rx');
                hold on
                plot(1:length(CnTemp(:)), CnTemp(:), '-k');
                legend('Used', 'Computed', 'interpreter', 'latex', ...
                    'FontSize', 14, 'NumColumns', 2);
                xlabel('$n$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$C_{n}(\beta_m)$', 'interpreter', ...
                    'latex', 'FontSize', 14);
                xlim([0, obj.niFD]);
                hold off 
            end               
        end
        function computeFDCm(obj, ShowPlot)
            % compute beta and eta values if not already populated
            if isempty(obj.etaFDStatic)
                if nargin > 1 && ShowPlot
                    computeEtaFD(obj, true);
                else
                    computeEtaFD(obj);
                end 
            end
            CmTemp = NaN*ones(1, obj.miFD);  
            for j = 1:obj.miFD
                CmTemp(j) = FourierCoefficientFD2(obj, obj.etaFDStatic(j));
            end
            cGet_ = abs(CmTemp) > obj.climFD;
            obj.FDCm = CmTemp(cGet_);
            obj.etaFD = obj.etaFDStatic(cGet_);
            obj.RFD = obj.RFDStatic(cGet_);
            cGetIdx = find(cGet_);                      
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(cGetIdx, obj.FDCm, 'rx');
                hold on
                plot(1:length(CmTemp(:)), CmTemp(:), '-k');
                legend('Used', 'Computed', 'interpreter', 'latex', ...
                    'FontSize', 14, 'NumColumns', 2);
                xlabel('$m$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$C_{m}(\eta_m)$', 'interpreter', ...
                    'latex', 'FontSize', 14);
                xlim([0, obj.miFD]);
                hold off 
            end               
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % top boundary velocity solution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeUbar(obj, ShowPlot)
            % computes velocity profile on top boundary
            obj.ubar = -obj.Qc./(2*pi*obj.rtop*obj.h);
            % plot velocity solutions on u-r axis to check solution
            obj.ubarfig = figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.4], 'Visible', 'off');
            plot(obj.rtop, obj.ubar, '.-k');
            title('Velocity Distribution in Top Boundary', ...
                'interpreter', 'latex', 'FontSize', 16);
            ylabel('$u/U_\infty$', 'interpreter', 'latex', 'FontSize', 16);
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 16);
            hold off 
            if nargin > 1 && ShowPlot
                 obj.ubarfig.Visible = 'on';
                 obj.ubarContFig.Visible = 'on';
            end 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % top boundary temperature distribution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function iterateThetaT(obj)
            % iterates the temperature for discretized energy equation
            n = length(obj.zhat); m = length(obj.rtop);
            % initialize thetaT if first time step
            if isempty(obj.FT0), obj.FT0 = ones(n, m); end
            % compute stagnant region temperature if not populated already
            if isempty(obj.thetaS), iterateThetaS(obj); end 
            % compute velocity distribution if not populated already
            if isempty(obj.ubar), computeUbar(obj); end 
            % compute matrix exponential if empty
            if isempty(obj.expTop), computeExpTop(obj); end
            % boundary condition at z=0 
            obj.thetaT(end, :) = obj.S2T;            
            % boundary condition at r = b
%             obj.thetaT(:, end) = obj.thetaS(ceil(length(obj.zbar)/2), end);
            obj.thetaT(:, end) = obj.thetaS(end, end);
            % set initial condition           
            thetaT_ = obj.thetaT';
            thetaT_ = obj.GT*thetaT_(:) + obj.LambdaT2*obj.thetaA;
            % calculate new temperature states with matrix exponential                    
            thetaT_ = obj.expTop*thetaT_(:);
            % set new temperatures for top region
            obj.thetaT = reshape(thetaT_, m, n)';
            % fix corner temperatures
            obj.thetaT(1, 1) = obj.thetaT(1, 2);
            obj.thetaT(end, 1) = obj.thetaT(end, 2);
            obj.thetaT(1, end) = obj.thetaT(2, end);
            obj.thetaT(end, end) = obj.thetaT(end-1, end);            
            % compute upper boundary temperature for center channel
            computeThetaChat(obj);
        end      
        function computeExpTop(obj)
            % computes the static top exponential
            n = length(obj.zhat); m = length(obj.rtop); nm = n*m;
            % top boundary meshed domain state matrix elements 
            obj.OmegaT1 = -2/obj.drtop^2 - 2/obj.dzhat^2;
            Ls = diag(ones(length(obj.ubar)-1, 1), -1);
            obj.OmegaT2 = 1/obj.drtop^2 + (1./repmat(obj.rtop', n, 1) ... 
                   - obj.Pe*repmat(Ls*obj.ubar', n, 1))/(2*obj.drtop);
            obj.OmegaT3 = 1/obj.dzhat^2;
            Us = diag(ones(length(obj.ubar)-1, 1), 1);
            obj.OmegaT4 = 1/obj.drtop^2 - (1./repmat(obj.rtop', n, 1) ... 
                   - obj.Pe*repmat(Us*obj.ubar', n, 1))/(2*obj.drtop);
            obj.OmegaT5 = 1/obj.dzhat^2;
            obj.AT = spdiags([obj.OmegaT1*ones(nm, 1), obj.OmegaT2, ...
                        obj.OmegaT3*ones(nm, 1), obj.OmegaT4, ...
                        obj.OmegaT5*ones(nm, 1)], [0, 1, m, -1, -m], ...
                        nm, nm);                                         
            % eliminate boundary elements from AT
            obj.AT(1:m, :) = 0;                 % z = h
            obj.AT(end-m+1:end, :) = 0;         % z = 0
            obj.AT(m:m:end, :) = 0;             % r = b
            obj.AT(1:m:end-m+1, :) = 0;         % r = a
            % Boundary state matrix
            obj.GT = eye(nm);
            % boundary condition at z=h
            obj.AT(1:m, :) = obj.AT(m+1:2*m, :)*1/(1 + obj.Bi5*obj.dzhat); 
            obj.LambdaT2 = [obj.Bi5*obj.dzhat/(1 + obj.Bi5*obj.dzhat) ...
                              *ones(m-1, 1); zeros(nm-m+1, 1)];
            obj.GT((m-1)*nm:nm+1:2*m*nm) = 1/(1 + obj.Bi5*obj.dzhat);
            obj.GT(1:nm+1:m*(nm+1)) = 0;
            % boundary condition at r = a
            obj.AT(1:m:end-m+1, :) = obj.AT(2:m:end-m+2, :);
            % matrix exponential for single time step
            obj.expTop = expm(obj.AT*obj.df);
        end
        function computeThetaChat(obj)
            % computes the homogeneous temperature for the top boundary in
            % the center channel based on an energy balance on the junction
            % that connects the top and center boundaries
            if isempty(obj.thetaT), iterateThetaT(obj);end
            if isempty(obj.ubar), computeUbar(obj); end
            if isempty(obj.wbar), computeWbar(obj); end
            fz = simpsonIntegrator(obj, obj.zhat);
            obj.thetaChat = 2/obj.a0*obj.ubar(1)/obj.wbar(1) ...
                                                    *fz*obj.thetaT(:, 1);           
        end
        function recordTopLoss(obj)
            % records dT from outer radius to inner radius in the top
            % channel at the current time step
            if isempty(obj.topLoss); obj.topLoss = cell(3, 1); end
            obj.topLoss{1} = [obj.topLoss{1}, obj.FoNow];
            obj.topLoss{2} = [obj.topLoss{2}, ...
                                  obj.thetaT(:, 2) - obj.thetaT(:, end-1)];
            obj.topLoss{3} = [obj.topLoss{3}, ...
                              obj.theta2T(obj.thetaT(:, 2)) ...
                            - obj.theta2T(obj.thetaT(:, end-1))];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center channel velocity distribution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function computeWbar(obj, ShowPlot)
            % compute bulk velocity distribution in center channel
            obj.wbar = -obj.Qc/(pi*(obj.a0)^2)*ones(size(obj.zcenter));
            % plot velocity solutions on u-z axis to check solutions
            obj.wbarfig = figure('Units', 'normalized', ...
                'Position', [0 0 0.3 0.6], 'Visible', 'off');
            plot(obj.wbar, obj.zcenter, '.-k');
            hold on
            title('Velocity Distribution in Center Channel', ...
                'interpreter', 'latex', 'FontSize', 16);
            xlabel('$w/U_\infty$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 16);
            ylim([0, max(obj.zcenter)]);
            hold off 
            if nargin > 1 && ShowPlot
                 obj.wbarfig.Visible = 'on';
            end 
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center boundary temperature distribution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function iterateThetaC(obj, sensitivity)
            % iterates the temperature for discretized energy equation
            n = length(obj.zcenter); m = length(obj.rhat); nm = n*m;
            if nargin < 2, sensitivity = 0; end
            % initialize thetaT if first time step
            if isempty(obj.FC0), obj.FC0 = ones(n, m); end
            % compute stagnant region temperature if not populated already
            if isempty(obj.thetaS), iterateThetaS(obj); end             
            % compute velocity distribution if not populated already
            computeWbar(obj);
            % compute top boundary temperature if not populated already
            if isempty(obj.thetaT), iterateThetaT(obj); end
            % top boundary meshed domain state matrix elements 
            obj.OmegaC1 = -2/obj.dzc^2 - 2/obj.drhat^2;
            obj.OmegaC2 = 1/obj.drhat^2 ...
                        + 1/(2*repmat(obj.rhat', n, 1)*obj.drhat);
            obj.OmegaC3 = 1/obj.dzc^2 + 0.1*obj.Pe ...
                   *repmat(obj.wbar', m, 1)/(2*obj.dzc);
            obj.OmegaC4 = 1/obj.drhat^2 ...
                        - 1/(2*repmat(obj.rhat', n, 1)*obj.drhat);
            obj.OmegaC5 = 1/obj.dzc^2 - 0.1*obj.Pe ...
                   *repmat(obj.wbar', m, 1)/(2*obj.dzc);
            obj.AC = spdiags([obj.OmegaC1*ones(nm, 1), obj.OmegaC2', ...
                              obj.OmegaC3, obj.OmegaC4', obj.OmegaC5], ...
                              [0, 1, m, -1, -m], nm, nm);
            % eliminate boundary elements from AC
            obj.AC(1:m, :) = 0;                 % z = 1
            obj.AC(end-m+1:end, :) = 0;         % z = 0
            obj.AC(m:m:end, :) = 0;             % r = a
            obj.AC(1:m:end-m+1, :) = 0;         % r = 0            
            % set boundary temperatures
            setCenterChannelBoundaries(obj);
            % compute new temperature states with matrix exponential
            thetaC_ = obj.thetaC';
            thetaC_ = expm(obj.AC*obj.df)*thetaC_(:); 
            % set new temperatures for top region
            obj.thetaC = reshape(thetaC_, m, n)';
            % compute sensitivity if wanted
            if sensitivity
                computeCenterSensitivity(obj);
            end
        end           
        function computeCenterSensitivity(obj)
            % computes sensativity w.r.t. select model parameters
            n = length(obj.zcenter); m = length(obj.rhat); nm = n*m;            
            if isempty(obj.scPe), obj.scPe = zeros(n, m); end
            if isempty(obj.scQc), obj.scQc = zeros(n, m); end
            if isempty(obj.sca), obj.sca = zeros(n, m); end
            % set initial conditions and 'constant' thetaC effects
            setdAC_Pe(obj);
            setdAC_Qc(obj);
            setdAC_a(obj);
            scPe0 = reshape(obj.scPe', [], 1);
            scQc0 = reshape(obj.scQc', [], 1);
            sca0 = reshape(obj.sca', [], 1);
            bPe = obj.dAC_Pe*reshape(obj.thetaC', [], 1);
            bQc = obj.dAC_Qc*reshape(obj.thetaC', [], 1);
            ba = obj.dAC_a*reshape(obj.thetaC', [], 1);
            % compute sensativity at current time step
            obj.scPe = expm(obj.df*[obj.AC, eye(nm); ...
                                    zeros(nm), zeros(nm)])*[scPe0; bPe];
            obj.scPe = reshape(abs(obj.scPe(1:nm)), m, n)';
            obj.scQc = expm(obj.df*[obj.AC, eye(nm); ...
                                    zeros(nm), zeros(nm)])*[scQc0; bQc];
            obj.scQc = reshape(abs(obj.scQc(1:nm)), m, n)';
            obj.sca = expm(obj.df*[obj.AC, eye(nm); ...
                                    zeros(nm), zeros(nm)])*[sca0; ba];
            obj.sca = reshape(abs(obj.sca(1:nm)), m, n)';                       
        end
        function setdAC_Pe(obj)
            % sets dAC_Pe matrix for sensitivity computation
            n = length(obj.zcenter); m = length(obj.rhat); nm = n*m;
            Omega3 = -repmat(obj.wbar', m, 1)/(2*obj.dzc);
            Omega5 = repmat(obj.wbar', m, 1)/(2*obj.dzc);
            obj.dAC_Pe = spdiags([Omega3, Omega5], [m, -m], nm, nm);
            % eliminate boundary elements from AT
            obj.dAC_Pe(1:m, :) = 0;                 % z = 1
            obj.dAC_Pe(end-m+1:end, :) = 0;         % z = 0
            obj.dAC_Pe(m:m:end, :) = 0;             % r = a
            obj.dAC_Pe(1:m:end-m+1, :) = 0;         % r = 0   
            % boundary condition at r = 0  
            obj.dAC_Pe(1:m:end-m+1, :) = obj.dAC_Pe(2:m:end-m+2, :);
            % boundary condition at z=0
            obj.dAC_Pe(end-m+1:end, :) = obj.dAC_Pe(end-2*m+1:end-m, :);            
        end
        function setdAC_Qc(obj)
            % sets dAC_Pe matrix for sensitivity computation
            n = length(obj.zcenter); m = length(obj.rhat); nm = n*m;
            Omega3 = obj.Pe/(2*obj.dzc*pi*obj.a0^2);
            Omega5 = -obj.Pe/(2*obj.dzc*pi*obj.a0^2);
            obj.dAC_Qc = spdiags([Omega3*ones(nm, 1), ...
                                  Omega5*ones(nm, 1)], [m, -m], nm, nm);
            % eliminate boundary elements from AT
            obj.dAC_Qc(1:m, :) = 0;                 % z = 1
            obj.dAC_Qc(end-m+1:end, :) = 0;         % z = 0
            obj.dAC_Qc(m:m:end, :) = 0;             % r = a
            obj.dAC_Qc(1:m:end-m+1, :) = 0;         % r = 0   
            % boundary condition at r = 0  
            obj.dAC_Qc(1:m:end-m+1, :) = obj.dAC_Qc(2:m:end-m+2, :);
            % boundary condition at z=0
            obj.dAC_Qc(end-m+1:end, :) = obj.dAC_Qc(end-2*m+1:end-m, :); 
        end
        function setdAC_a(obj)
            % sets dAC_Pe matrix for sensitivity computation
            n = length(obj.zcenter); m = length(obj.rhat); nm = n*m;
            Omega3 = -obj.Pe*obj.Qc/(obj.dzc*pi*obj.a0^3);
            Omega5 = obj.Pe*obj.Qc/(obj.dzc*pi*obj.a0^3);
            obj.dAC_a = spdiags([Omega3*ones(nm, 1), ...
                                  Omega5*ones(nm, 1)], [m, -m], nm, nm);
            % eliminate boundary elements from AT
            obj.dAC_a(1:m, :) = 0;                 % z = 1
            obj.dAC_a(end-m+1:end, :) = 0;         % z = 0
            obj.dAC_a(m:m:end, :) = 0;             % r = a
            obj.dAC_a(1:m:end-m+1, :) = 0;         % r = 0   
            % boundary condition at r = 0  
            obj.dAC_a(1:m:end-m+1, :) = obj.dAC_a(2:m:end-m+2, :);
            % boundary condition at z=0
            obj.dAC_a(end-m+1:end, :) = obj.dAC_a(end-2*m+1:end-m, :); 
        end
        function setCenterChannelBoundaries(obj)
            % sets center boundary conditions
            m = length(obj.rhat);
            % boundary condition at r = 0  
            obj.AC(1:m:end-m+1, :) = obj.AC(2:m:end-m+2, :);
            obj.thetaC(:, 1) = obj.thetaC(:, 2);
            % boundary condition at z=0
            obj.AC(end-m+1:end, :) = obj.AC(end-2*m+1:end-m, :);
            obj.thetaC(end, 2:end) = obj.thetaC(end-1, 2:end);
            % boundary condition at r = a
            obj.thetaC(:, end) = obj.S2C; 
            % boundary condition at z=1
            obj.thetaC(1, :) = obj.thetaChat;
        end               
        function theta_ = computeThetaO(obj)
            % computes the centerline and bulk outlet temperature for 
            % current time
%             if isempty(obj.thetaC)
%                 iterateThetaC(obj);
%             end
%             if isempty(obj.thetaO)
%                 obj.thetaO = [];
%             end
%             obj.thetaO = [obj.thetaO; [obj.FoNow, obj.thetaC(end, 1)]];
%             obj.thetaOB = [obj.thetaOB; [obj.FoNow, ...
%                 mean(obj.thetaC(end, :))]];
            obj.thetaO = []; obj.thetaOB = [];
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii), 1);
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                r_ = obj.theta{ii_, 3};
                [~, ri] = min(abs(r_ - obj.a0));
                fr = simpsonIntegrator(obj, r_(1:ri));
                theta_  = obj.theta{ii_, 1}(1:ri);
                obj.thetaO = [obj.thetaO; [t, theta_(1)]];
                obj.thetaOB = [obj.thetaOB; [t, ...
                                   2/r_(ri)^2*(theta_.*r_(1:ri))*fr']];
                               
            end
            theta_ = obj.thetaOB;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % storage bin quiescent air temperature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeThetaA(obj, thetaP, z_, r_)
            % computes the air temperature at the current iteration
            t = obj.Fo2t(obj.df, 1); n = length(z_); m = length(r_);
            % heat loss (W) at the particle boundary
            if isempty(obj.qTopP)
                obj.qTopP = cell(2, 1);
                obj.qTopP{1} = zeros(1, m);
            else
                obj.qTopP{1} = obj.qTopP{2};
            end
            fr = simpsonIntegrator(obj, r_);
            obj.qTopP{2} = -obj.kp*(thetaP(end, :) - thetaP(end-1, :))./ ...
                         (z_(end) - z_(end-1))*(obj.T0 - obj.Tinf)/obj.Hp;
            qTopP_ = 2*pi*fr*(obj.qTopP{1}.*r_)'*obj.Hp^2;
            qTopP__ = 2*pi*fr*((obj.qTopP{2} - obj.qTopP{1}) ...
                                                     ./t.*r_)'*obj.Hp^2; 
            % heat loss (W) at the top bin boundary
            At = pi*obj.bp^2*obj.Hp^2;
            qTop_ = obj.qLossT(1)*At;
            qTop__ = (obj.qLossT(2) - obj.qLossT(1))/t*At;
            % heat loss (W) at wall boundary
            [~, zwi] = min(abs(obj.ztop - obj.zbarW));
            zw = obj.zbarW(zwi:end);
            fz = simpsonIntegrator(obj, zw);
            qWall_ = 2*pi*obj.bp*fz*obj.qWall{1, 1}(zwi:end)*obj.Hp^2;
            qWall__ = 2*pi*obj.bp*fz*(obj.qWall{2, 1}(zwi:end) - ...
                                   obj.qWall{1, 1}(zwi:end))/t*obj.Hp^2;
            % compute thetaA
            cp_ = 1005;     % (J/kgK)
            rho_ = 1.225;   % (m3/kg)
            V = pi*obj.b^2*(1 - obj.ztop)*obj.Hp^3;
            C_ = cp_*rho_*V;            
            TA = obj.theta2T(obj.thetaA) + t*(qTopP_ - qTop_ - qWall_)/C_ ...
                + 0.5*t^2*(qTopP__ - qTop__ - qWall__)/C_;
            obj.thetaA = obj.T2theta(TA);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % storage bin top, base and wall RC model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeUbase(obj)
            % computes the overall heat transfer coefficient for the base
            % using the insulation cell array
            R = 0;
            for i = 1:size(obj.baseInsulation, 1)
                t = abs(obj.baseInsulation{i, 2}(2) - ...
                                            obj.baseInsulation{i, 2}(1));
                k_ = obj.baseInsulation{i, 3};
                R = R + t/(pi*(obj.bp*obj.Hp)^2*k_);
            end
            R = R + 1/(pi*(obj.bp*obj.Hp)^2*obj.hInf);            
            obj.hp2 = 1/(pi*(obj.bp*obj.Hp)^2*R);
        end
        function computeUwall(obj)
            % computes the overall heat transfer coefficient for the wall
            % using the insulation cell array
            R = 0;
            for i = 1:size(obj.wallInsulation, 1)
                r1 = obj.wallInsulation{i, 2}(1);
                r2 = obj.wallInsulation{i, 2}(2);
                k_ = obj.wallInsulation{i, 3};
                R = R + log(r2/r1)/(2*pi*obj.Hp*k_);
            end
            R = R + 1/(2*pi*r2*obj.Hp*obj.hInf);
            obj.hp4 = 1/(2*pi*obj.bp*obj.Hp^2*R);
        end
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
            % initialize state variables and boundary condition
            preheatTime = 3600*1;
            y0 = computeBaseSys(obj, [1, 1], [0, obj.t2Fo(preheatTime, 1)]);
            obj.thetaBase = y0(end, N+1:end);
%             computeBaseSys(obj);
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
                    obj.Cwall(i, i+1) = -(obj.T0 - obj.Tinf)/ ...
                        ((obj.Rwall{i, 2} + obj.Rcw{i, 2}) ...
                        *2*pi*obj.meshedWallInsulation{i, 2}(1)*obj.Hp);
                end
            end
            obj.Dwall = zeros(2*N, 1);
            obj.wallSys = ...
                ss(full(obj.Awall), obj.Bwall, obj.Cwall, obj.Dwall);
            % initialize state variables and boundary condition (with preheat) 
            preheatTime = 3600*10;
            y0 = computeWallSys(obj, [1, 1], [0, obj.t2Fo(preheatTime, 1)]);
            obj.thetaWall = y0(end, N+1:end);
%             computeWallSys(obj);
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
            N = size(obj.meshedBaseInsulation, 1);
            y = lsim(obj.baseSys, u, Fo_, obj.thetaBase);
            obj.qLossB = y(end, 1); obj.thetaBase = y(end, N+1:end);
            obj.g2 = -obj.qLossB*obj.Hp/(obj.kp*(obj.T0 - obj.Tinf));
%             obj.g2 = obj.g2*(1 - exp(-obj.FoNow./obj.tauW1));
        end
        function y = computeWallSys(obj, u, Fo_)
            % computes the heat flux leaving the wall
            if nargin < 2, u = ones(2, 1); end
            if nargin < 3, Fo_ = linspace(0, obj.df, length(u)); end
            N = size(obj.meshedWallInsulation, 1);
            y = lsim(obj.wallSys, u, Fo_, obj.thetaWall);
            obj.qLossW = y(end, 1); obj.thetaWall = y(end, N+1:end);
            obj.g4 = obj.qLossW*obj.Hp/(obj.kp*(obj.T0 - obj.Tinf));
%             obj.g4 = obj.g4*(1 - exp(-obj.FoNow./obj.tauW1));
        end
        function y = computeWallSensitivityR(obj, u, yss, Fo_)
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
        function setLumpedBaseMesh(obj)
            % adds discrete lumped elements for the wall model to increase
            % spatial resolution
            N = size(obj.baseInsulation, 1);
            obj.meshedBaseInsulation = {};
            p_ = 1;
            for i = 1:N
                M = obj.baseInsulation{i, 6};
                zmin = obj.baseInsulation{i, 2}(1);
                zmax = obj.baseInsulation{i, 2}(2);
                z_ = linspace(zmin, zmax, M+1);
                for j = 1:M
                    obj.meshedBaseInsulation{p_, 1} = ...
                        obj.baseInsulation{i, 1};
                    obj.meshedBaseInsulation{p_, 2} = [z_(j), z_(j+1)];
                    obj.meshedBaseInsulation{p_, 3} = ...
                        obj.baseInsulation{i, 3};
                    obj.meshedBaseInsulation{p_, 4} = ...
                        obj.baseInsulation{i, 4};
                    obj.meshedBaseInsulation{p_, 5} = ...
                        obj.baseInsulation{i, 5};
                    if j == 1
                        obj.meshedBaseInsulation{p_, 6} = ...
                            obj.baseInsulation{i, 7};
                    else
                        obj.meshedBaseInsulation{p_, 6} = 1e16;
                    end
                    p_ = p_+1;
                end             
            end 
            obj.zBL = unique(cell2mat(obj.meshedBaseInsulation(:, 2)), 'sorted');
        end
        function setLumpedWallMesh(obj)
            % adds discrete lumped elements for the wall model to increase
            % spatial resolution
            N = size(obj.wallInsulation, 1);
            obj.meshedWallInsulation = {};
            p_ = 1;
            for i = 1:N
                M = obj.wallInsulation{i, 6};
                rmin = obj.wallInsulation{i, 2}(1);
                rmax = obj.wallInsulation{i, 2}(2);
                r_ = linspace(rmin, rmax, M+1);
                for j = 1:M
                    obj.meshedWallInsulation{p_, 1} = ...
                        obj.wallInsulation{i, 1};
                    obj.meshedWallInsulation{p_, 2} = [r_(j), r_(j+1)];
                    obj.meshedWallInsulation{p_, 3} = ...
                        obj.wallInsulation{i, 3};
                    obj.meshedWallInsulation{p_, 4} = ...
                        obj.wallInsulation{i, 4};
                    obj.meshedWallInsulation{p_, 5} = ...
                        obj.wallInsulation{i, 5};
                    if j == 1
                        obj.meshedWallInsulation{p_, 6} = ...
                        obj.wallInsulation{i, 7};
                    else
                        obj.meshedWallInsulation{p_, 6} = 1e16;
                    end
                    p_ = p_+1;
                end             
            end 
            obj.rWL = unique(cell2mat(obj.meshedWallInsulation(:, 2)), 'sorted');
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
        function dae = erDAE(obj, t, y)
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
        % storage bin base and wall 1D continuous model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        function computeThetaW1D(obj, Fo_)
            % computes the bin wall temperature distribution at time Fo_
            % using the stored (or default) initial condition and domain
            % knowledge
            if isempty(obj.thetaW), buildWall(obj); end
            if nargin < 2, Fo_ = obj.Fo(end); end         
            M = size(obj.wallInsulation, 1);
            % initial condition contribution
            computeThetaWIC1D(obj, Fo_);
            % boundary condition contribution
            if isempty(obj.gW1), obj.gW1 = obj.Biw1*1; end
            obj.gW1t = [obj.gW1t, obj.gW1(1)];
            [~, k_] = min(abs(obj.Fo - Fo_));
            computeThetaWBC1D(obj, obj.Fo(1:k_));
            % Green's formula
            for i = 1:M
                obj.thetaW{i} = obj.thetaWIC{i} + obj.thetaWBC{i};
            end
            % compute heat rates through wall
            ki = [obj.wallInsulation{:, 3}];
            obj.qWall = -ki(M)*(obj.theta2T(obj.thetaW{M}(:, end)) ...
                      - obj.theta2T(obj.thetaW{M}(:, end-1)))/ ...
                      (obj.Hp*(obj.rbarW{M}(end) - obj.rbarW{M}(end - 1)));
            obj.qLossW = -ki(1)*(obj.theta2T(obj.thetaW{1}(:, 2)) ...
                      - obj.theta2T(obj.thetaW{1}(:, 1)))/ ...
                      (obj.Hp*(obj.rbarW{M}(2) - obj.rbarW{M}(1)));                         
        end  
        function computeThetaWIC1D(obj, Fo_)
            % computes the initial condition contribution for the Green's
            % formula for the temperature distribution through the wall          
            if isempty(obj.GWCm), computeGWCm(obj); end
            M = size(obj.wallInsulation, 1);
            n = length(obj.zbarW); m = cell(M, 1); 
            for i = 1:M, m{i} = length(obj.rbarW{i}); end
            k_ = length(Fo_);             
            % compute Fourier series for filtered solution
            for i = 1:M, obj.thetaWIC{i} = sparse(zeros(n*k_, m{i}*k_)); end
            % compute temperature solution for each composite layer
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            for i = 1:M
                for j = 1:length(obj.GWCnm)
                    bWm_ = obj.bWm{j};
                    A = bWm_(2*i-1); B = bWm_(2*i);
                    C = obj.GWCnm(j);
                    F = XWt(obj, Fo_, 0, obj.etaW(j), alphai(i));
                    Z = ones(n, 1); 
                    R = XWm(obj, obj.rbarW{i}, obj.etaW(j), A, B);                                                    
                    obj.thetaWIC{i} = obj.thetaWIC{i} + ...
                                           sparse(kron(C*diag(F), Z'*R));                   
                end
            end
        end
        function computeThetaWBC1D(obj, Fo_, ShowPlot)
            % computes the inner-wall boundary condition contribution for
            % the Green's formula for the temperature distribution through
            % the wall layers    
            % compute eta values if not already populated
            if isempty(obj.etaWStatic)
                if nargin > 2 && ShowPlot
                    computeEtaW(obj, true);   
                else
                    computeEtaW(obj);
                end     
            end
            M = size(obj.wallInsulation, 1);
            n = length(obj.zbarW); m = cell(M, 1);
            for i = 1:M
                m{i} = length(obj.rbarW{i}); 
                obj.thetaWBC{i} = zeros(n, m{i});
            end
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            ki = [obj.wallInsulation{:, 3}];
            for i = 1:M
                if obj.miWBCtot > length(obj.etaWStatic) 
                    mtot = length(obj.etaWStatic);
                else
                    mtot = obj.miWBCtot;
                end
                for m = 1:mtot
                    eta_ = obj.etaWStatic(m);
                    [A_, c] = computeAWm(obj, eta_);
                    b_ = A_(1:end-1, 2:end)\c;
                    bWm_ = [1; b_];
                    Ai = bWm_(2*i-1); Bi = bWm_(2*i); % need to make sure that the constants for each composite layer are coming from the same j (as this is what the boundary conditions are satisfied with)
                    C = ki(1)/alphai(1)*XWm(obj, obj.rbarW{1}(1), ...
                        eta_, bWm_(1), bWm_(2))/NWm(obj, eta_, bWm_); 
                    if length(Fo_) == 1
                        C = 0;
                    elseif length(Fo_) == 2
                        C = trapz(Fo_, C*obj.gW1t.*XWt(obj, -Fo_, 0, ...
                                                         eta_, alphai(i)));
                    else
                        ft = simpsonIntegrator(obj, Fo_);
                        C = (C*obj.gW1t.*XWt(obj, -Fo_, 0, eta_, ...
                                                        alphai(i)))*ft';
                    end
                    if C < 1e-5, C = 0; end
                    F = XWt(obj, Fo_(end), 0, eta_, alphai(i));
                    Zi = ones(n, 1);
                    Ri = XWm(obj, obj.rbarW{i}, eta_, Ai, Bi);                                                    
                    obj.thetaWBC{i} = obj.thetaWBC{i} + C*F*Zi*Ri;  
                end
            end                     
        end
        function c = FourierCoefficientGW1D(obj, eta_, bWm_)
            % computes fourier coefficient for given inputs for the
            % initial condition eigenvalue problem of the composite wall
            Cnum = computeThetaWICIntegral1D(obj, eta_, bWm_);
            c = Cnum/NWm(obj, eta_, bWm_);
        end
        function x = computeThetaWICIntegral1D(obj, eta_, bWm_)
            % computes integral from numerator of Fourier coefficient for
            % the filtered thetaW solution
            M = size(obj.wallInsulation, 1);
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            ki = [obj.wallInsulation{:, 3}];
            x = 0;
            for i = 1:M
                r_ = obj.rbarW{i};
                A = bWm_(2*i-1); B = bWm_(2*i);
                alpha_ = alphai(i); rhoW_ = obj.rhoW{i};
                fr = simpsonIntegrator(obj, r_);
                Ir = (XWm(obj, r_, eta_, A, B).*r_.*rhoW_(1, :))*fr';
                x = x + ki(i)/alpha_*Ir;
            end                       
        end
        function computeGWCm(obj, ShowPlot)
            % compute new beta values
            if isempty(obj.etaWStatic)
                if nargin > 1 && ShowPlot
                    computeEtaW(obj, true);
                else
                    computeEtaW(obj);
                end
            end        
            M = size(obj.wallInsulation, 1);
            if obj.miW > length(obj.etaWStatic) 
                miW_ = length(obj.etaWStatic); 
            else
                miW_ = obj.miW;
            end
            for n = 1:M
                CnmTemp = NaN*ones(1, miW_);
                etaTemp = NaN*ones(1, miW_); 
                bWmTemp = cell(1, miW_);
                for j = 1:miW_
                    CnmTemp(j) = FourierCoefficientGW1D(obj, ...
                        obj.etaWStatic(j), obj.bWmStatic{j});
                    etaTemp(j) = obj.etaWStatic(j);
                    bWmTemp(j) = obj.bWmStatic(j);
                end
                obj.cGet = abs(CnmTemp) > obj.climW;
                Cnm_ = CnmTemp(obj.cGet);
    %             w = obj.betaW*pi/max(obj.beta); 
    %             sigmaN = sin(w)./w;
                sigmaN = ones(length(Cnm_), 1);
                obj.GWCnm = sigmaN.*Cnm_;
                obj.etaW = etaTemp(obj.cGet);
                obj.bWm = bWmTemp(obj.cGet);
                cGetIdx = find(obj.cGet);  
                % plot to check solutions
                if nargin > 1 && ShowPlot
                    figure('Units', 'normalized', ...
                        'Position', [0 0 0.4 0.25], 'Visible', 'on');
                    plot(cGetIdx, Cnm_, 'rx');
                    hold on
                    plot(1:length(CnmTemp(:)), CnmTemp(:), '-k');
                    legend('Used', 'Computed', 'interpreter', 'latex', ...
                        'FontSize', 14, 'NumColumns', 2);
                    xlabel('$nm$', 'interpreter', 'latex', 'FontSize', 14);
                    ylabel('$C_{m}(\eta_m)$', 'interpreter', ...
                        'latex', 'FontSize', 14);
                    xlim([0, miW_]);
                    hold off 
                end    
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % storage bin base and wall 2D continuous model (old)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeThetaW(obj, Fo_) 
            % computes the bin wall temperature distribution at time Fo_
            % using the stored (or default) initial condition and domain
            % knowledge
            if isempty(obj.thetaW), buildWall(obj); end
            if nargin < 2, Fo_ = obj.Fo(end); end         
            M = size(obj.wallInsulation, 1);
            % initial condition contribution
            computeThetaWIC(obj, Fo_);
            % boundary condition contribution
            if isempty(obj.gW1), obj.gW1 = obj.Biw1*ones(size(obj.zbarW)); end
%             obj.gW1t = [obj.gW1t, obj.gW1];
            [~, k_] = min(abs(obj.Fo - Fo_));
            computeThetaWBC(obj, obj.Fo(1:k_));
            % Green's formula
            for i = 1:M
                obj.thetaW{i} = obj.thetaWIC{i} + obj.thetaWBC{i};
            end
            % compute heat rates through wall
            ki = [obj.wallInsulation{:, 3}];
            obj.qWall = -ki(M)*(obj.theta2T(obj.thetaW{M}(:, end)) ...
                      - obj.theta2T(obj.thetaW{M}(:, end-1)))/ ...
                      (obj.Hp*(obj.rbarW{M}(end) - obj.rbarW{M}(end - 1)));
            obj.qLossW = -ki(1)*(obj.theta2T(obj.thetaW{1}(:, 2)) ...
                      - obj.theta2T(obj.thetaW{1}(:, 1)))/ ...
                      (obj.Hp*(obj.rbarW{M}(2) - obj.rbarW{M}(1)));                         
        end        
        function computeThetaWTemp(obj, Fo_)
            % computes the temperature distribution through the wall layers           
            if isempty(obj.thetaW), buildWall(obj); 
            else, obj.rhoW = obj.thetaW; end
            M = size(obj.wallInsulation, 1);
            % compute steady-state filter and coefficients
            computeFW(obj); computeCnmW(obj);
            if nargin < 2, Fo_ = obj.Fo; end
            n = length(obj.zbarW); m = cell(M, 1); 
            for i = 1:M, m{i} = length(obj.rbarW{i}); end
            k_ = length(Fo_);             
            % compute Fourier series for filtered solution
            for i = 1:M, obj.thetaW{i} = sparse(zeros(n*k_, m{i}*k_)); end
            % compute temperature solution for each composite layer
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            ki = [obj.wallInsulation{:, 3}];
            if obj.FoModeNow == 'H', z_ = obj.zbarH; 
            else, z_ = obj.zbar; end
            obj.g4 = zeros(length(z_), 1);
            n = length(obj.zbarW);
            if isempty(obj.qWall)
                obj.qWall = cell(2, M+1);
%                 obj.qWall{1, 1} = zeros(length(obj.zbarW), 1);
                obj.qWall{1, 1} = sparse(zeros(n*k_, m{1}*k_));
            else
                obj.qWall{1, 1} = obj.qWall{2, 1};
            end
%             obj.qWall{2, 1} = zeros(length(obj.zbarW), 1);
            obj.qWall{2, 1} = sparse(zeros(n*k_, k_));
            for i = 1:M
%                 m = length(obj.rbarW{i}); 
                if isempty(obj.qWall{1, i+1})
%                     obj.qWall{1, i+1} = zeros(length(obj.zbarW), 1);
                    obj.qWall{1, i+1} = sparse(zeros(n*k_, k_));
                else
                    obj.qWall{1, i+1} = obj.qWall{2, i+1};
                end
%                 obj.qWall{2, i+1} = zeros(length(obj.zbarW), 1);
                obj.qWall{2, i+1} = sparse(zeros(n*k_, k_));
                bWm_ = obj.bWm{i};
                % compute for current time step
%                 obj.thetaW{i} = zeros(n, m);
                for j = 1:length(obj.CnmW{i})
                    A = bWm_{j}(2*i-1); B = bWm_{j}(2*i);
                    C = obj.CnmW{i}(j);
                    F = Xt(obj, Fo_, obj.betaW{i}(j), obj.etaW{i}(j));
                    Z = XWn(obj, obj.zbarW, obj.betaW{i}(j), alphai(i)); 
                    R = XWm(obj, obj.rbarW{i}, obj.etaW{i}(j), ...
                                                        A, B, alphai(i));                                                    
%                     obj.thetaW{i} = obj.thetaW{i} + C*F*Z'*R;
                    obj.thetaW{i} = obj.thetaW{i} + ...
                                           sparse(kron(C*diag(F), Z'*R));
                    % compute heat loss at outer radius of section                     
                    [~, Rp] = XWm(obj, obj.rbarW{i}(end), ...
                                    obj.etaW{i}(j), A, B, alphai(i));
%                     obj.qWall{2, i+1} = obj.qWall{2, i+1} + ...
%                                   ki(i)*(obj.T0 - obj.Tinf)/obj.Hp* ...
%                                   C*F*Z'*Rp;
                    obj.qWall{2, i+1} = obj.qWall{2, i+1} + ...
                                  ki(i)*(obj.T0 - obj.Tinf)/obj.Hp* ...
                                  sparse(kron(C*diag(F), Z'*Rp));    
                    if i == 1                       
                        [~, Rp] = XWm(obj, obj.rbarW{1}(1), ...
                                        obj.etaW{1}(j), A, B, alphai(1));
%                         obj.qWall{2, 1} = obj.qWall{2, 1} + ...
%                                       ki(1)*(obj.T0 - obj.Tinf)/obj.Hp* ...
%                                       C*F*Z'*Rp; 
                        obj.qWall{2, 1} = obj.qWall{2, 1} + ...
                                      ki(1)*(obj.T0 - obj.Tinf)/obj.Hp* ...
                                      sparse(kron(C*diag(F), Z'*Rp)); 
                        Z = XWn(obj, z_, obj.betaW{i}(j), alphai(i)); 
                        obj.g4 = obj.g4 + ki(1)/obj.kp*C*F(end)*Z'*Rp;
                    end                    
                end
                obj.thetaW{i} = obj.thetaW{i} + ...
                                   sparse(kron(eye(length(Fo_)), obj.FW{i}));
                obj.rhoW{i} = obj.thetaW{i}((k_-1)*n+1:end, (k_-1)*m{i}+1:end);
                if i == 1
                    obj.g4 = -(obj.g4 + obj.g4f); 
%                     obj.qWall{2, 1} = -(obj.qWall{2, 1} + obj.qWallf{1});
                    obj.qWall{2, 1} = -(obj.qWall{2, 1} + ...
                                         sparse(kron(eye(length(Fo_)), ...
                                         obj.qWallf{1})));
                end
%                 obj.qWall{2, i+1} = -(obj.qWall{2, i+1} + obj.qWallf{i+1});
                obj.qWall{2, i+1} = -(obj.qWall{2, i+1} + ...
                                         sparse(kron(eye(length(Fo_)), ...
                                         obj.qWallf{i+1})));
            end
        end
        function computeThetaWIC(obj, Fo_)
            % computes the initial condition contribution for the Green's
            % formula for the temperature distribution through the wall          
%             if isempty(obj.GWCnm), computeGWCnm(obj); end
            M = size(obj.wallInsulation, 1);
            n = length(obj.zbarW); m = cell(M, 1); 
            for i = 1:M, m{i} = length(obj.rbarW{i}); end
            k_ = length(Fo_);             
            % compute Fourier series for filtered solution
            for i = 1:M, obj.thetaWIC{i} = sparse(zeros(n*k_, m{i}*k_)); end
            % compute temperature solution for each composite layer
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            for i = 1:M
                for j = 1:length(obj.GWCnm)
                    bWm_ = obj.bWm{j};
                    A = bWm_(2*i-1); B = bWm_(2*i);
                    C = obj.GWCnm(j);
                    F = XWt(obj, Fo_, obj.betaW(j), obj.etaW(j), alphai(i));
                    Z = XWn(obj, obj.zbarW, obj.betaW(j)); 
                    R = XWm(obj, obj.rbarW{i}, obj.etaW(j), A, B);                                                    
                    obj.thetaWIC{i} = obj.thetaWIC{i} + ...
                                           sparse(kron(C*diag(F), Z'*R));                   
                end
            end
        end
        function computeThetaWBC(obj, Fo_, ShowPlot)
            % computes the inner-wall boundary condition contribution for
            % the Green's formula for the temperature distribution through
            % the wall layers    
            % compute new beta values
            if isempty(obj.betaWStatic)
                if nargin > 2 && ShowPlot
                    computeBetaW(obj, true);
                else
                    computeBetaW(obj);
                end
            end
            % compute eta values if not already populated
            if isempty(obj.etaWStatic)
                if nargin > 2 && ShowPlot
                    computeEtaW(obj, true);
                else
                    computeEtaW(obj);
                end     
            end
            M = size(obj.wallInsulation, 1);
            n = length(obj.zbarW); m = cell(M, 1);
            for i = 1:M
                m{i} = length(obj.rbarW{i}); 
                obj.thetaWBC{i} = zeros(n, m{i});
            end
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            ki = [obj.wallInsulation{:, 3}];
%             fz = simpsonIntegrator(obj, obj.zbarW);
            for i = 1:M
                if obj.miWBCtot > length(obj.etaWStatic) 
                    mtot = length(obj.etaWStatic);
                else
                    mtot = obj.miWBCtot;
                end
                if obj.niWBCtot > length(obj.betaWStatic) 
                    ntot = length(obj.betaWStatic);
                else
                    ntot = obj.niWBCtot;
                end
                for m = 1:mtot
                    eta_ = obj.etaWStatic(m);
                    [A_, c] = computeAWm(obj, eta_);
                    b_ = A_(1:end-1, 2:end)\c;
                    bWm_ = [1; b_];
                    for n = 1:ntot
                        beta_ = obj.betaWStatic(n);
                        Ai = bWm_(2*i-1); Bi = bWm_(2*i); % need to make sure that the constants for each composite layer are coming from the same j (as this is what the boundary conditions are satisfied with)
                        C = ki(1)/alphai(1)*XWm(obj, obj.rbarW{1}(1), ...
                            eta_, bWm_(1), bWm_(2)) ...
                            /(NWn(obj, beta_)*NWm(obj, eta_, bWm_)); 
%                         Zj = XWn(obj, obj.zbarW, beta_);
%                         IBC = fz*(C*obj.gW1.*Zj');
                        if m == 1
                            Iz = computeZWintNow(obj, beta_);
                            obj.IBCW{n} = [obj.IBCW{n}, C*Iz];
                        end
                        if length(Fo_) == 1
                            C = 0;
                        elseif length(Fo_) == 2
                            C = trapz(Fo_, obj.IBCW{n}.*XWt(obj, -Fo_', beta_, ...
                                                  eta_, alphai(i)));
                        else
                            ft = simpsonIntegrator(obj, Fo_);
                            C = (obj.IBCW{n}.*XWt(obj, -Fo_', beta_, eta_, ...
                                                    alphai(i)))*ft';
                        end
                        if C < 1e-5, C = 0; end
                        F = XWt(obj, Fo_(end), beta_, eta_, alphai(i));
                        Zi = XWn(obj, obj.zbarW, beta_);
                        Ri = XWm(obj, obj.rbarW{i}, eta_, Ai, Bi);                                                    
                        obj.thetaWBC{i} = obj.thetaWBC{i} + C*F*Zi'*Ri;  
                    end
                end
            end                     
        end
        function Iz = computeZWintNow(obj, beta_)
            % computes the analytical integral for the inner boundary
            % Green's function at the current time step by treating the
            % wall function as a step function (particle region and air
            % region)
            [~, zi] = min(abs(obj.ztop - obj.zbarW));
            if obj.ztop == 0, g1_ = 0; else, g1_ = obj.gW1(1:zi); end
            if obj.ztop == 1, g2_ = 0; else, g2_ = obj.gW1(zi+1:end); end             
            if beta_ < 1e-15
                Iz = mean(g1_)*obj.ztop + mean(g2_)*(1 - obj.ztop);
            else
                Iz = mean(g1_)/beta_*sin(beta_*obj.ztop) + ...
                     mean(g2_)/beta_*(sin(beta_) - sin(beta_*obj.ztop));
            end
        end
        function computeThetaWBCOld(obj, Fo_)
            % computes the inner-wall boundary condition contribution for
            % the Green's formula for the temperature distribution through
            % the wall layers            
            computeGWBCCnm(obj, 1);
            M = size(obj.wallInsulation, 1);
            n = length(obj.zbarW); m = cell(M, 1);
            for i = 1:M
                m{i} = length(obj.rbarW{i}); 
                obj.thetaWBC{i} = zeros(n, m{i});
            end
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);          
            for i = 1:M
                for j = 1:length(obj.GWBCCnm)
                    bWm_ = obj.bWBCm{j};
                    A = bWm_{j}(2*i-1); B = bWm_{j}(2*i); % need to make sure that the constants for each composite layer are coming from the same j (as this is what the boundary conditions are satisfied with)
                    C = obj.GWBCCnm(j); 
                    F = Xt(obj, Fo_(end), obj.betaWBC(j), ...
                                                obj.etaWBC(j), alphai(i));
                    Z = XWn(obj, obj.zbarW, obj.betaWBC(j));
                    R = XWm(obj, obj.rbarW{i}, obj.etaWBC(j), A, B);
                    obj.thetaWBC{i} = obj.thetaWBC{i} + C*F*Z'*R;                 
                end
            end                     
        end      
        function c = FourierCoefficientGW(obj, beta_, eta_, bWm_)
            % computes fourier coefficient for given inputs for the
            % initial condition eigenvalue problem of the composite wall
            Cnum = computeThetaWICIntegral(obj, beta_, eta_, bWm_);
            c = Cnum/(NWn(obj, beta_)*NWm(obj, eta_, bWm_));
        end
        function c = FourierCoefficientGWBC(obj, beta_, eta_, bWm_)
            % computes the fourier coefficient for given inputs for the 
            % boundary value problem of the composite wall
            Cnum = computeThetaWBCIntegral(obj, beta_, eta_, bWm_);
            M = size(obj.wallInsulation, 1); c = cell(1, M);
            for j = 1:M
                c{j} = Cnum{j}/(NWn(obj, beta_)*NWm(obj, eta_, bWm_));
            end
        end
        function x = computeThetaWICIntegral(obj, beta_, eta_, bWm_)
            % computes integral from numerator of Fourier coefficient for
            % the filtered thetaW solution
            M = size(obj.wallInsulation, 1);
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            ki = [obj.wallInsulation{:, 3}];
            x = 0;
            z_ = obj.zbarW;
            for i = 1:M
                r_ = obj.rbarW{i};
                A = bWm_(2*i-1); B = bWm_(2*i);
                alpha_ = alphai(i); rhoW_ = obj.rhoW{i};
                fr = simpsonIntegrator(obj, r_);
                fz = simpsonIntegrator(obj, z_);
                Ir = (XWm(obj, r_, eta_, A, B).*r_.*rhoW_ )*fr';
                x = x + ...
                    ki(i)/alpha_*(XWn(obj, z_, beta_).*Ir')*fz';
            end                       
        end
        function x = computeThetaWBCIntegral(obj, beta_, eta_, bWm_)
            % computes integral from numerator of Fourier coefficient for
            % the filtered thetaW solution            
            M = size(obj.wallInsulation, 1);
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            ki = [obj.wallInsulation{:, 3}];
            x = cell(1, M);
            for j = 1:M
                fz = simpsonIntegrator(obj, obj.zbarW);
                C = ki(j)/alphai(j)*XWm(obj, ...
                       obj.rbarW{1}(1), eta_, bWm_(1), bWm_(2));
                Z = XWn(obj, obj.zbarW, beta_);
                x{j} = fz*(C*obj.gW1t.*Z');
%                 if length(Fo_) == 1
%                     x = 0;
%                 elseif length(Fo_) == 2
%                     x = x + trapz(Fo_, IBC.*XWt(obj, -Fo_, beta_, ...
%                                                         eta_, alphai(i)));
%                 else
%                     ft = simpsonIntegrator(obj, Fo_);
%                     x = x + (IBC.*XWt(obj, -Fo_, beta_, eta_, ...
%                                                           alphai(i)))*ft';
%                 end  
            end
        end       
        function computeGWCnm_old(obj, ShowPlot)
            % compute new beta values
            if isempty(obj.betaWStatic)
                if nargin > 1 && ShowPlot
                    computeBetaW(obj, true);
                else
                    computeBetaW(obj);
                end
            end
            % compute eta values if not already populated
            if isempty(obj.etaWStatic)
                if nargin > 1 && ShowPlot
                    computeEtaW(obj, true);
                else
                    computeEtaW(obj);
                end     
            end           
            M = size(obj.wallInsulation, 1);
            if obj.miW > length(obj.etaWStatic) 
                    miW_ = length(obj.etaWStatic); 
            else
                    miW_ = obj.miW;
            end
            for n = 1:M
                if obj.niW > length(obj.betaWStatic) 
                    niW_ = length(obj.betaWStatic); 
                else
                    niW_ = obj.niW;
                end
                CnmTemp = NaN*ones(niW_, miW_);
                betaTemp = NaN*ones(niW_, miW_);
                etaTemp = NaN*ones(niW_, miW_); 
                bWmTemp = cell(niW_, miW_);
                for i = 1:niW_
                    for j = 1:miW_
                        CnmTemp(i, j) = FourierCoefficientGW(obj, ...
                            obj.betaWStatic(i), obj.etaWStatic(j), ...
                            obj.bWmStatic{j});
                        betaTemp(i, j) = obj.betaWStatic(i);
                        etaTemp(i, j) = obj.etaWStatic(j);
                        bWmTemp(i, j) = obj.bWmStatic(j);
                    end
                end
                obj.cGet = abs(CnmTemp) > obj.climW;
                Cnm_ = CnmTemp(obj.cGet);
                obj.betaW = betaTemp(obj.cGet);
    %             w = obj.betaW*pi/max(obj.beta); 
    %             sigmaN = sin(w)./w;
                sigmaN = ones(length(obj.betaW), 1);
                obj.GWCnm = sigmaN.*Cnm_;
                obj.etaW = etaTemp(obj.cGet);
                obj.bWm = bWmTemp(obj.cGet);
                cGetIdx = find(obj.cGet);  
                % plot to check solutions
                if nargin > 1 && ShowPlot
                    figure('Units', 'normalized', ...
                        'Position', [0 0 0.4 0.25], 'Visible', 'on');
                    plot(cGetIdx, Cnm_, 'rx');
                    hold on
                    plot(1:length(CnmTemp(:)), CnmTemp(:), '-k');
                    legend('Used', 'Computed', 'interpreter', 'latex', ...
                        'FontSize', 14, 'NumColumns', 2);
                    xlabel('$nm$', 'interpreter', 'latex', 'FontSize', 14);
                    ylabel('$C_{nm}(\beta_n, \eta_m)$', 'interpreter', ...
                        'latex', 'FontSize', 14);
                    xlim([0, niW_*miW_]);
                    hold off 
                end    
            end
        end
        function computeGWBCCnm(obj, ShowPlot)
            % compute new beta values
            if isempty(obj.betaWStatic)
                if nargin > 2 && ShowPlot
                    computeBetaW(obj, true);
                else
                    computeBetaW(obj);
                end
            end
            % compute eta values if not already populated
            if isempty(obj.etaWStatic)
                if nargin > 2 && ShowPlot
                    computeEtaW(obj, true);
                else
                    computeEtaW(obj);
                end     
            end           
            if obj.miWBC > length(obj.etaWStatic) 
                    miW_ = length(obj.etaWStatic); 
            else
                    miW_ = obj.miWBC;
            end
            if obj.niWBC > length(obj.betaWStatic) 
                niW_ = length(obj.betaWStatic); 
                obj.niWBC = niW_;
            else
                niW_ = obj.niWBC;
            end  
            M = size(obj.wallInsulation, 1);
            CnmTemp = cell(niW_, miW_);
            CnmMag = zeros(niW_, miW_);
            betaTemp = NaN*ones(niW_, miW_);
            etaTemp = NaN*ones(niW_, miW_); 
            bWmTemp = cell(niW_, miW_);
            for i = 1:niW_
                for j = 1:miW_
                    CnmTemp{i, j} = FourierCoefficientGWBC(obj, ...
                        obj.betaWStatic(i), obj.etaWStatic(j), ...
                        obj.bWmStatic{j});
                    betaTemp(i, j) = obj.betaWStatic(i);
                    etaTemp(i, j) = obj.etaWStatic(j);
                    bWmTemp(i, j) = obj.bWmStatic(j);
                    for n = 1:M 
                        CnmMag(i, j) = CnmMag(i, j) + CnmTemp{i, j}{n};
                    end
                end
            end
            obj.cGet = abs(CnmMag) > M*obj.climWBC;
            Cnm_ = CnmTemp(obj.cGet);
            obj.betaWBC = betaTemp(obj.cGet);
%             w = obj.betaW*pi/max(obj.beta); 
%             sigmaN = sin(w)./w;
%             sigmaN = ones(length(obj.betaWBC), 1);
            obj.GWBCCnm = Cnm_;
            obj.etaWBC = etaTemp(obj.cGet);
            obj.bWBCm = bWmTemp(obj.cGet);
            cGetIdx = find(obj.cGet);  
            % plot to check solutions
            if nargin > 1 && ShowPlot
                for n = 1:M
                    figure('Units', 'normalized', ...
                        'Position', [0 0 0.4 0.25], 'Visible', 'on');
                    Cget = zeros(length(Cnm_), 1);
                    for i = 1:length(Cnm_)
                        Cget(i) = Cnm_{i}{n};                        
                    end
                    plot(cGetIdx, Cget, 'rx'); hold on;
                    idx = 1:length(CnmTemp(:));
                    C_ = CnmTemp(:);
                    Cplot = zeros(length(CnmTemp(:)), 1);
                    for i = 1:length(CnmTemp(:))
                            Cplot(i) = C_{i}{n}; 
                    end
                    plot(idx, Cplot, '-k')
                    legend('Used', 'Computed', 'interpreter', 'latex', ...
                        'FontSize', 14, 'NumColumns', 2);
                    xlabel('$nm$', 'interpreter', 'latex', 'FontSize', 14);
                    ylabel('$CBC_{nm}(\beta _n, \eta _m)$', 'interpreter', ...
                        'latex', 'FontSize', 14);
                    xlim([0, niW_*miW_]);
                    hold off 
                end
            end    
        end       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % storage bin base and wall 2D continuous model (new)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeThetaW_new(obj, Fo_)
            % uses a Green's function heat kernel to compute the current
            % 2D wall temperature
            if isempty(obj.thetaW), buildWall(obj); end
            if isempty(obj.betaW), computeBetaW(obj); end
            if isempty(obj.etaW), computeEtaW(obj); end
%             if isempty(obj.GWCnm), computeGWCnmSimilarity(obj); end
            if isempty(obj.gW1)
                obj.gW1 = obj.Biw1*ones(length(obj.zbarW), 1); 
            end
            fz = simpsonIntegrator(obj, obj.zbarW);
            M = size(obj.wallInsulation, 1);
            k1 = obj.wallInsulation{1, 3};
            alpha1 = [obj.wallInsulation{1, 3}]./ ...
                ([obj.wallInsulation{1, 4}].*[obj.wallInsulation{1, 5}]);
            % compute wall heat kernel for current time
            if isempty(obj.KWI) 
                obj.KWI = cell(M, length(obj.betaW), length(obj.etaW)); 
                obj.KW = cell(M, length(obj.betaW), length(obj.etaW)); 
            end
            for i = 1:M
                alphai = [obj.wallInsulation{i, 3}]./ ...
                ([obj.wallInsulation{i, 4}].*[obj.wallInsulation{i, 5}]);                
                obj.thetaW{i} = ...
                    zeros(length(obj.zbarW), length(obj.rbarW{i}));
                for n = 1:length(obj.betaW)
                    for m = 1:length(obj.etaW)
                        bWm_ = obj.bWm{m};
                        Ai = bWm_(2*i-1); Bi = bWm_(2*i);
                        A1 = bWm_(1); B1 = bWm_(2);
                        beta_ = obj.betaW(n); eta_ = obj.etaW(m);
                        lambda_ = sqrt(beta_^2 + eta_^2);
                        if Fo_ == 0
                            obj.KWI{i, n, m} = 0; obj.KW{i, n, m} = 0;
                        else
                            if beta_ > 0
                                obj.KWI{i, n, m} = [obj.KWI{i, n, m}, ...
                                    (exp(alphai/obj.alphapPacked ...
                                    *lambda_^2*Fo_) ...
                                    *XWn(obj, obj.zbarW, beta_) ...
                                    .*obj.gW1')*fz'];
                            else
                                obj.KWI{i, n, m} = [obj.KWI{i, n, m}, 0];
                            end
                            if length(obj.KWI{i, n, m}) == 2
                                obj.KW{i, n, m} = trapz(obj.Fo(1:2), obj.KWI{i, n, m});
                            else
                                ft = simpsonIntegrator(obj, ...
                                    obj.Fo(1:length(obj.KWI{i, n, m})));
                                obj.KW{i, n, m} = obj.KWI{i, n, m}*ft';
                            end
                        end
                        F = exp(-alphai/obj.alphapPacked*lambda_^2*Fo_);
                        C = XWm(obj, obj.rbarW{1}(1), eta_, A1, B1) ...
                            *obj.KW{i, n, m}/(obj.betaWN(n)*obj.etaWN(m));
                        R = XWm(obj, obj.rbarW{i}, eta_, Ai, Bi);
                        Z = XWn(obj, obj.zbarW, beta_);
                        obj.thetaW{i} = obj.thetaW{i} + C*F*Z'*R;
                    end
                end
                obj.thetaW{i} = alphai*k1/(obj.alphapPacked*alpha1) ...
                                                        *obj.thetaW{i};
            end            
        end        
        function xi = XWn(~, z_, beta_)
            % computes z-BVP eigenfunctions for all composite layers with
            % the corresponding properties (alpha_)
            xi = cos(beta_*z_);
        end
        function [xi, xip] = XWm(~, r_, eta_, A, B)
            % computes r-BVP eigenfunctions for all composite layers with
            % the corresponding coefficients defined in A and B
            if eta_ < 1e-15
                xi = A + B*log(r_); xip = B./r_;
            else
                xi = A*besselj(0, eta_*r_) + B*bessely(0, eta_*r_);
                xip = -eta_*(A*besselj(1, eta_*r_) + B*bessely(1, eta_*r_));
            end
        end
        function t = XWt(obj, Fo_, beta_, eta_, alphai)
            % transient component of analytic solutions
            t = exp(-alphai/obj.alphaPacked*(beta_^2 + eta_^2)*Fo_);
        end
        function [x, y, xp, yp] = phiW(~, r_, eta_)
            % r-dimension eigenfunctions with eigenvalue, eta_, and thermal
            % diffusivity alpha_
            if eta_ < 1e-15
                x = 1; y = log(r_); xp = 0; yp = 1/r_;
            else
                x = besselj(0, eta_*r_); y = bessely(0, eta_*r_);
                xp = -eta_.*besselj(1, eta_*r_); 
                yp = -eta_.*bessely(1, eta_*r_);
            end
        end  
        function ni = NWn(obj, beta_)
            % z-BVP eigenvalue problem norm
            M = size(obj.wallInsulation, 1);
            if beta_ < 1e-14, ni = M; else, ni = 0.5*M; end
        end
        function ni = NWm(obj, eta_, bWm_)
            % r-BVP eigenvalue problem norm
            M = size(obj.wallInsulation, 1);
            ri = unique([obj.wallInsulation{:, 2}])/obj.Hp;
            alphai = [obj.wallInsulation{:, 3}]./ ...
                ([obj.wallInsulation{:, 4}].*[obj.wallInsulation{:, 5}]);
            ki = [obj.wallInsulation{:, 3}];
            ni = 0;
            for i = 1:M
                A = bWm_(2*i-1); B = bWm_(2*i);
                r1 = ri(i); r2 = ri(i+1); 
                q1 = r1*eta_; 
                q2 = q1*r2/r1;
                ni = ni + ki(i)/alphai(i)*(0.5*A^2* ...
                    (r2^2*(besselj(0, q2)^2 + besselj(1, q2)^2) ...
                    - r1^2*(besselj(0, q1)^2 + besselj(1, q1)^2)) + ...
                    0.5*B^2*(r2^2*(bessely(0, q2)^2 + bessely(1, q2)^2) ...
                    - r1^2*(bessely(0, q1)^2 + bessely(1, q1)^2)) + ...
                    A*B*(r2^2*(besselj(0, q2)*bessely(0, q2) + ...
                    besselj(1, q2)*bessely(1, q2)) - ...
                    r1^2*(besselj(0, q1)*bessely(0, q1) + ...
                    besselj(1, q1)*bessely(1, q1))));
            end
        end   
        function zi = ZWn(~, beta_)
            % transcendental equation for finding r-dimension eigenvalues
            zi = sin(beta_*1);                       
        end
        function zi = ZWm(obj, eta_)
            % transcendental equation for finding r-dimension eigenvalues
            [A, ~] = computeAWm(obj, eta_); zi = det(A);
%             if rcond(A) < 1e-15, zi = 0; end
            if ~isreal(zi), zi = 1; end         
        end
        function [A, c] = computeAWm(obj, eta_)
            % computes the r-BVP boundary condition matrix
            M = size(obj.wallInsulation, 1);
            A = zeros(2*M);
            c = zeros(2*M-1, 1);
            ri = unique([obj.wallInsulation{:, 2}])/obj.Hp;
            ki = [obj.wallInsulation{:, 3}];
            % fill boundary at inner surface
            [x, y, xp, yp] = phiW(obj, ri(1), eta_);
            A(1, 1:2) = [-xp + obj.Biw1*x, -yp + obj.Biw1*y];
            c(1) = -(-xp + obj.Biw1*x);
            % fill boundaries at composite connections
            for i = 1:M-1
                [x1, y1, xp1, yp1] = phiW(obj, ri(i+1), eta_);
                [x2, y2, xp2, yp2] = phiW(obj, ri(i+1), eta_);
                A(2*i:2*i+1, 2*i-1:2*i+2) = [x1, y1, -x2, -y2; ...
                     ki(i)*xp1, ki(i)*yp1, -ki(i+1)*xp2, -ki(i+1)*yp2]; 
                if i == 1, c(2:3) = [-x1; -ki(i)*xp1]; end 
            end
            % fill exterior surface boundary
            [x, y, xp, yp] = phiW(obj, ri(end), eta_);
            A(end, end-1:end) = [xp + obj.Biw2*x, yp + obj.Biw2*y];                        
        end
        function computeBetaW(obj, ShowPlot)
            interval = linspace(0, obj.pfW, obj.pW);   % interval/spacing 
                                                     % of root calculation
            rn = NaN*ones(obj.pW, 1);              % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.pW
                rn(i) = fzero(@(beta_) ...
                    ZWn(obj, beta_), interval(i), options);
            end
            obj.betaW = rn(diff(rn)>1e-10); % only keep unique roots
            obj.betaW(1) = 0;
            obj.betaWN = obj.betaW;
            for i = 1:length(obj.betaW) 
                obj.betaWN(i) = NWn(obj, obj.betaW(i));
            end       
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', 'Position', [0 0 0.4 0.25]);
                plot(interval, ZWn(obj, interval), '-k');
                hold on
                plot(obj.betaW, ZWn(obj, obj.betaW), '.r');
                legend('$Z_n(\hat{\beta})$', '$\beta_n$', 'interpreter', ...
                    'latex', 'FontSize', 14, 'NumColumns', 2);
                xlabel('$\hat{\beta}$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylabel('$Z_n(\hat{\beta})$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                xlim([0, obj.pfW]);
                hold off 
            end                   
        end                  
        function computeEtaW(obj, ShowPlot)
            % computes r-BVP eigenvalues, corresponding boundary
            % condition matrices, and eigenfunction coefficients
            interval = linspace(1, obj.qfW, obj.qW);   % interval/spacing 
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
            obj.etaWN = obj.etaW;
            for i = 1:length(obj.etaW)                            
                [A_, c] = computeAWm(obj, obj.etaW(i));
                b_ = A_(2:end, 2:end)\c;
                obj.bWm{i} = [1; b_]; 
                obj.etaWN(i) = NWm(obj, obj.etaW(i), obj.bWm{i});
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
        function computeGWCnmSimilarity(obj, ShowPlot)
            % compute new beta values
            if isempty(obj.betaW)
                if nargin > 1 && ShowPlot
                    computeBetaW(obj, true);
                else
                    computeBetaW(obj);
                end
            end
            % compute eta values if not already populated
            if isempty(obj.etaW)
                if nargin > 1 && ShowPlot
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
            if obj.niW > length(obj.betaW) 
                niW_ = length(obj.betaW); 
            else
                niW_ = obj.niW;
            end
            k1 = obj.wallInsulation{1, 3};
            alpha1 = [obj.wallInsulation{1, 3}]./ ...
                ([obj.wallInsulation{1, 4}].*[obj.wallInsulation{1, 5}]);
            CnmTemp = NaN*ones(niW_, miW_);
            betaTemp = NaN*ones(niW_, miW_);
            etaTemp = NaN*ones(niW_, miW_); 
            bWmTemp = cell(niW_, miW_);
            for i = 1:niW_
                for j = 1:miW_
                    if obj.betaW(i) == 0
                        CnmTemp(i, j) = 0;
                    else
                        CnmTemp(i, j) = k1/alpha1*XWm(obj, obj.rbarW{1}(1), ...
                            obj.etaW(j), obj.bWm{j}(1), obj.bWm{j}(2))/ ...
                            (obj.betaWN(i)*obj.etaWN(j));
                    end
                    betaTemp(i, j) = obj.betaW(i);
                    etaTemp(i, j) = obj.etaW(j);
                    bWmTemp(i, j) = obj.bWm(j);
                end
            end
            obj.cGet = abs(CnmTemp) > obj.climW;
            Cnm_ = CnmTemp(obj.cGet);
            obj.betaW = betaTemp(obj.cGet);
%             w = obj.betaW*pi/max(obj.beta); 
%             sigmaN = sin(w)./w;
            sigmaN = ones(length(obj.betaW), 1);
            obj.GWCnm = sigmaN.*Cnm_;
            obj.etaW = etaTemp(obj.cGet);
            obj.bWm = bWmTemp(obj.cGet);
            cGetIdx = find(obj.cGet);  
            % plot to check solutions
            if nargin > 1 && ShowPlot
                figure('Units', 'normalized', ...
                    'Position', [0 0 0.4 0.25], 'Visible', 'on');
                plot(cGetIdx, Cnm_, 'rx');
                hold on
                plot(1:length(CnmTemp(:)), CnmTemp(:), '-k');
                legend('Used', 'Computed', 'interpreter', 'latex', ...
                    'FontSize', 14, 'NumColumns', 2);
                xlabel('$nm$', 'interpreter', 'latex', 'FontSize', 14);
                ylabel('$C_{nm}(\beta_n, \eta_m)$', 'interpreter', ...
                    'latex', 'FontSize', 14);
                xlim([0, niW_*miW_]);
                hold off 
            end    
        end
        function buildWall(obj)
            % assembles wall parameters from the materials defined in
            % wallInsulation
            obj.zbarW = nodeGen(obj, [0, 1], obj.nzbarW);
            for i = 1:size(obj.wallInsulation, 1)
                obj.rbarW{i} = nodeGen(obj, ...
                    obj.wallInsulation{i, 2}./obj.Hp, obj.nrbarW{i});
                obj.thetaW{i} = ...
                    zeros(length(obj.zbarW), length(obj.rbarW{i}));
            end
            obj.Biw1 = obj.hcw*obj.Hp/obj.wallInsulation{1, 3};
            obj.Biw1A = obj.hcwA*obj.Hp/obj.wallInsulation{1, 3};
            obj.Biw2 = obj.hInf*obj.Hp/obj.wallInsulation{end, 3};
            if isempty(obj.rhoW), obj.rhoW = obj.thetaW; end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Domain modification or data filling
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function updateDomain(obj)
            % updates domain according to mass accounting
            storeTotalMass(obj);
            obj.HNow = obj.HNow - obj.dH;
            obj.ztop = obj.HNow/obj.H;
            fitzMesh(obj);
            fitTheta(obj);
        end
        function updateDomainC(obj)
            % updates domain according to mass accounting
            storeTotalMass(obj);
            obj.HNow = obj.HNow + obj.dHCh;
            obj.ztop = obj.HNow/obj.H;
            fitzMeshC(obj);
            fitThetaC(obj);
        end
        function fitzMeshC(obj)
            % fits z meshes to be between 0 and ztop, resizing dzbar and
            % dzc accordingly. Adjusts the length according to the time
            % step relative to the set update frequency, modZ.
%             [~, n] = min(abs(obj.FoNow - obj.Fo));
%             if mod(n, obj.modZH) == 0
%                 obj.nzH = length(obj.zbarH) + 1; 
%             else
%                 obj.nzH = length(obj.zbarH);
%             end
%             [obj.zbarH, obj.dzH] = ...
%                               nodeGen(obj, [0, obj.ztop], obj.nzH);  
            obj.nzH = ceil(obj.nzH0*obj.ztop);
            [obj.zbarH, obj.dzH] = ...
                              nodeGen(obj, [0, obj.ztop], obj.nzH);
        end
        function fitzMesh(obj)
            % fits z meshes to be between 0 and ztop, resizing dzbar and
            % dzc accordingly. Adjusts the length according to the time
            % step relative to the set update frequency, modZ.
%             [~, n] = min(abs(obj.FoNow - obj.Fo));
%             if mod(n, obj.modZS) == 0
%                 obj.nzbar = length(obj.zbar) - 1; 
%             else
%                 obj.nzbar = length(obj.zbar);
%             end
%             if mod(n, obj.modZC) == 0
%                 nzc = length(obj.zcenter) - 1;
%             else
%                 nzc = length(obj.zcenter);
%             end
%             [obj.zbar, obj.dzbar] = ...
%                             nodeGen(obj, [0, obj.ztop], obj.nzbar);
%             obj.zcenter = linspace(0, obj.ztop, nzc);
%             obj.dzc = obj.zcenter(2);
            obj.nzbar = ceil(obj.nzbar0*obj.ztop);
            obj.nzc = ceil(obj.nzc0*obj.ztop);
            [obj.zbar, obj.dzbar] = ...
                            nodeGen(obj, [0, obj.ztop], obj.nzbar);
            [obj.zcenter, obj.dzc] = nodeGen(obj, [0, obj.ztop], obj.nzc);
            computeWbar(obj);
            computeUbar(obj);            
        end
        function fitThetaC(obj)
            % fits/updates theta matrices to match current mesh
            move = abs(size(obj.thetaH, 1) - length(obj.zbarH));
            if move > 0             
                % update temperature 
                for n = 1:move
                    obj.thetaH = [obj.thetaH; ...
                           obj.thetaCi*ones(1, length(obj.thetaH(1, :)))];
                    % compute temperature distribution for mixing depth
%                     setMixingProfileC(obj);
                    % set new initial condition for next time step
                    obj.rhoH = obj.thetaH; 
                end
            end         
        end
        function setMixingProfileC(obj)
            % computes bulk temperature/s for the mixing depth for the 
            % current time step in charging mode
            [~, im] = min(abs(obj.zbarH - (obj.ztop - obj.deltaM)));
            zm = obj.zbarH(im:end);
            if length(zm) < 3
                im = im - 3;
                zm = obj.zbarH(im:end);
            end
            obj.thetaMC = mean(obj.thetaH(im:end, :));
            obj.thetaH(im:end, :) = repmat(obj.thetaMC, length(zm), 1);
        end
        function fitTheta(obj)
            % fits/updates theta matrices to match current mesh
            move = abs(size(obj.thetaS, 1) - length(obj.zbar));
            movec = abs(size(obj.thetaC, 1) - length(obj.zcenter));
            if move > 0             
                % update temperature 
                for n = 1:move
                    for i = 1:length(obj.zhat)-1
                        obj.thetaT(i, :) = obj.thetaT(i+1, :);                    
                    end
                    obj.thetaT(end, :) = obj.S2T;
                    obj.thetaS = obj.thetaS(1:end-1, :); 
                    obj.FS0 = obj.FS0(1:end-1, :);
                end
            end
            if movec > 0
                for n = 1:movec
                    obj.thetaC = obj.thetaC(1:end-1, :);
                    obj.scPe = obj.scPe(1:end-1, :);
                    obj.scQc = obj.scQc(1:end-1, :);
                    obj.sca = obj.sca(1:end-1, :);                    
                end
            end             
        end
        function patchTheta(obj, k_, verify_conservation, sensitivity)
            % combines thetaS, thetaT and thetaC for time step k_
            if nargin < 3, verify_conservation = 0; end
            if nargin < 4, sensitivity = 0; end
            ns = length(obj.zhat); nl = length(obj.zbar); n = ns + nl - 1;
            ms = length(obj.rhat); ml = length(obj.rbar); m = ms + ml - 1;
            theta_ = NaN*ones(n, m);
            % save and reinitialize theta if no more space
            if mod(k_-1, obj.ls) == 0 && k_ ~= 1
                saveTheta(obj, k_-1);
                obj.theta = {};                
            end
            % normalize time step
            kn_ = mod(k_-1, obj.ls) + 1;
            % overlay current temperature in stagnant region
            theta_(1:nl, ms:m) = obj.thetaS;
            % overlay current temperature in top region
            thetaTMatch = obj.T2S('full');
            theta_(nl+1:n, ms:m) = flipud(thetaTMatch(1:end-1, :));
            % overlay current temperature in center channel
            thetaCMatch = obj.C2S('full');
            theta_(1:nl, 1:ms-1) = flipud(thetaCMatch(:, 1:end-1));
            % overlay bulk junction temperature in top-center junction
            theta_(nl+1:n, 1:ms-1) = mean(obj.thetaChat);
            % compute heat loss at wall, base, and top
            
            % append to theta cell array that will be saved to drive
            obj.theta{kn_, 1} = theta_;
            % store corresponding mesh arrays
            obj.theta{kn_, 2} = [obj.zbar, ...
                                 obj.zbar(end) + obj.zhat(2:end)];
            obj.theta{kn_, 3} = [obj.rhat(1:end-1), obj.rbar]; 
            obj.theta{kn_, 4} = k_;
            obj.theta{kn_, 5} = obj.Fo(k_);
            obj.theta{kn_, 6} = obj.FoMode{1, k_};
            obj.theta{kn_, 7} = obj.thetaWall;
            obj.theta{kn_, 8} = obj.thetaBase;
            obj.theta{kn_, 9} = obj.thetaRoof;
            obj.theta{kn_, 10} = obj.thetaWA;
            obj.theta{kn_, 11} = obj.qWall;
            obj.theta{kn_, 12} = obj.thetaA;
            [qw, qb, qt, qTot] = computeBinHeatLoss(obj, ...
                theta_(:, end-1:end), ...
                obj.theta{kn_, 3}(end) - obj.theta{kn_, 3}(end-1), obj.theta{kn_, 2}, ...
                theta_(1:2, :), ...
                obj.theta{kn_, 2}(2) - obj.theta{kn_, 2}(1), obj.theta{kn_, 3}, ...
                theta_(end-2:end-1, :), ...
                obj.theta{kn_, 2}(end-1) - obj.theta{kn_, 2}(end-2), obj.theta{kn_, 3});
            obj.theta{kn_, 13} = qw; %obj.qLossW;
            obj.theta{kn_, 14} = qb; %obj.qLossB;
            obj.theta{kn_, 15} = qt; %obj.qLossT;
            obj.theta{kn_, 16} = qTot; %obj.qTopP;
            obj.theta{kn_, 17} = obj.qRadR;
            obj.theta{kn_, 18} = obj.qRadW;
            obj.theta{kn_, 19} = obj.T2thetaK(obj.Tp);
            obj.theta{kn_, 20} = 0;
            if verify_conservation
                computeEnergyD(obj);
                obj.theta{kn_, 21} = obj.energy;
            end
            if sensitivity
                obj.theta{kn_, 22} = obj.scPe;
                obj.theta{kn_, 23} = obj.scQc;
                obj.theta{kn_, 24} = obj.sca;
            end
            % save if last time step
            if k_ == length(obj.Fo)
                saveTheta(obj, k_);
            end
        end
        function patchThetaC(obj, k_, verify_conservation)
            % save thetaH array and domain info
            if mod(k_-1, obj.ls) == 0 && k_ ~= 1
                saveTheta(obj, k_-1);
                obj.theta = {};                
            end
            % normalize time step
            kn_ = mod(k_-1, obj.ls) + 1;                        
            % append to theta cell array that will be saved to drive
            obj.theta{kn_, 1} = obj.thetaH;
            % compute heat loss at wall, base, and top
            [qw, qb, qt, qTot] = computeBinHeatLoss(obj, ...
                obj.thetaH(:, end-1:end), ...
                obj.rbarH(end) - obj.rbarH(end-1), obj.zbarH, ...
                obj.thetaH(1:2, :), ...
                obj.zbarH(2) - obj.zbarH(1), obj.rbarH, ...
                obj.thetaH(end-2:end-1, :), ...
                obj.zbarH(end-1) - obj.zbarH(end-2), obj.rbarH);
            % store corresponding mesh arrays
            obj.theta{kn_, 2} = obj.zbarH;
            obj.theta{kn_, 3} = obj.rbarH; 
            obj.theta{kn_, 4} = k_;
            obj.theta{kn_, 5} = obj.Fo(k_);
            obj.theta{kn_, 6} = obj.FoMode{1, k_};
            obj.theta{kn_, 7} = obj.thetaWall;
            obj.theta{kn_, 8} = obj.thetaBase;
            obj.theta{kn_, 9} = obj.thetaRoof;
            obj.theta{kn_, 10} = obj.thetaWA;
            obj.theta{kn_, 11} = obj.qWall;
            obj.theta{kn_, 12} = obj.thetaA;
            obj.theta{kn_, 13} = qw; %obj.qLossW;
            obj.theta{kn_, 14} = qb; %obj.qLossB;
            obj.theta{kn_, 15} = qt; %obj.qLossT;
            obj.theta{kn_, 16} = qTot; %obj.qTopP;
            obj.theta{kn_, 17} = obj.qRadR;
            obj.theta{kn_, 18} = obj.qRadW;
            obj.theta{kn_, 19} = obj.T2thetaK(obj.Tp);
            obj.theta{kn_,20} = obj.QChp*obj.rhopPack*obj.Fo2t(obj.df, 1)*obj.cpp*(obj.T0 - obj.Tinf)/1000;
            if verify_conservation
                computeEnergyH(obj);
                obj.theta{kn_, 21} = obj.energy;
            end
            % save if last time step
            if k_ == length(obj.Fo)
                saveTheta(obj, k_);
            end
        end
        function patchThetaH(obj, k_, verify_conservation)
            % save thetaH array and domain info
            if mod(k_-1, obj.ls) == 0 && k_ ~= 1
                saveTheta(obj, k_-1);
                obj.theta = {};                
            end
            % normalize time step
            kn_ = mod(k_-1, obj.ls) + 1;                        
            % append to theta cell array that will be saved to drive
            obj.theta{kn_, 1} = obj.thetaH;
            % compute heat loss at wall, base, and top
            [qw, qb, qt, qTot] = computeBinHeatLoss(obj, ...
                obj.thetaH(:, end-1:end), ...
                obj.rbarH(end) - obj.rbarH(end-1), obj.zbarH, ...
                obj.thetaH(1:2, :), ...
                obj.zbarH(2) - obj.zbarH(1), obj.rbarH, ...
                obj.thetaH(end-2:end-1, :), ...
                obj.zbarH(end-1) - obj.zbarH(end-2), obj.rbarH);
            % store corresponding mesh arrays
            obj.theta{kn_, 2} = obj.zbarH;
            obj.theta{kn_, 3} = obj.rbarH; 
            obj.theta{kn_, 4} = k_;
            obj.theta{kn_, 5} = obj.Fo(k_);
            obj.theta{kn_, 6} = obj.FoMode{1, k_};
            obj.theta{kn_, 7} = obj.thetaWall;
            obj.theta{kn_, 8} = obj.thetaBase;
            obj.theta{kn_, 9} = obj.thetaRoof;
            obj.theta{kn_, 10} = obj.thetaWA;
            obj.theta{kn_, 11} = obj.qWall;
            obj.theta{kn_, 12} = obj.thetaA;
            obj.theta{kn_, 13} = qw; %obj.qLossW;
            obj.theta{kn_, 14} = qb; %obj.qLossB;
            obj.theta{kn_, 15} = qt; %obj.qLossT;
            obj.theta{kn_, 16} = qTot; %obj.qTopP;
            obj.theta{kn_, 17} = obj.qRadR;
            obj.theta{kn_, 18} = obj.qRadW;
            obj.theta{kn_, 19} = obj.T2thetaK(obj.Tp);
            obj.theta{kn_, 20} = 0;
            if verify_conservation
                computeEnergyH(obj);
                obj.theta{kn_, 21} = obj.energy;
            end
            % save if last time step
            if k_ == length(obj.Fo)
                saveTheta(obj, k_);
            end
        end 
        function patchThetaK(obj, k_, thetaP)
            % compiles and saves data cell that contains results for the
            % simulation that matches Kevin's model 
            rW_ = []; thetaW_ = []; M = length(obj.rbarW);
            n = length(obj.zbarW); m = cell(M, 1); thetaWC_ = cell(M, 1);
            for i = 1:M, m{i} = length(obj.rbarW{i}); end
            for i = 1:length(obj.thetaW)
                rW_ = [rW_, obj.rbarW{i}];
%                 thetaWC_{i} = obj.thetaW{i}((k_-1)*n+1:k_*n, ...
%                                                  (k_-1)*m{i}+1:k_*m{i});
                thetaW_ = [thetaW_, obj.thetaW{i}];
%                 thetaW_ = [thetaW_, thetaWC_{i}];
            end
            n = length(obj.zbarW);
            mp = length(obj.rbarH); mw = length(rW_); m = mp + mw;
            theta_ = NaN*ones(n, m);
            % save and reinitialize theta if no more space
            if mod(k_-1, obj.ls) == 0 && k_ ~= 1
                saveThetaK(obj, k_-1);
                obj.theta = {};                
            end
            % normalize time step
            kn_ = mod(k_-1, obj.ls) + 1;
            % include particle/air temperature
            theta_(1:n, 1:mp) = thetaP;
            % include wall temperature
            theta_(1:n, mp+1:end) = thetaW_;
            % store data
            obj.thetaK{kn_, 1} = theta_;
            obj.thetaK{kn_, 2} = obj.zbarW;
            obj.thetaK{kn_, 3} = [obj.rbarH, rW_]; 
            obj.thetaK{kn_, 4} = k_;
            obj.thetaK{kn_, 5} = obj.Fo(k_);
            obj.thetaK{kn_, 6} = obj.thetaW;
            obj.thetaK{kn_, 7} = obj.qWall;
            obj.thetaK{kn_, 8} = obj.thetaA;
            obj.thetaK{kn_, 9} = obj.qLossW;
            % save if last time step
            if k_ == length(obj.Fo)
                saveThetaK(obj, k_);
            end
        end
        function computeThetaI(obj)
            % computes the average temperature of the top and center flow
            % channel boundaries to set the offset temperature that avoids
            % boundary osscilations in the Green's function
            tr = bulkTempR(obj, mean(obj.thetaT, 1), obj.rtop);
            tz = bulkTempZ(obj, mean(obj.thetaC, 2));
            obj.thetaI = mean([tr, tz]);
        end
        function savePoints(obj)
            % saves temperature value at discrete points defined in zStore
            % and rStore
            if isempty(obj.thetaStore)                
                for i = 1:length(obj.zStore)
                    c = struct('z', obj.zStore(i), 'r', obj.rStore(i), ...
                               'Fo', obj.Fo, 'theta', []);
                    obj.thetaStore{i} = c;
                end
            end
            for i = 1:length(obj.zStore)
                obj.thetaStore{i}.theta = [obj.thetaStore{i}.theta, ...
                    obj.thetaZR(obj.zStore(i), obj.rStore(i))];
            end
        end
        function x = thetaZR(obj, z, r)
            % extracts current nondimensional temperature at (z, r) from
            % the thetaS, thetaT and thetaC arrays
            if r < obj.a0
                % pull from thetaC
                [~, i] = min(abs(obj.zcenter - z));
                [~, j] = min(abs(obj.rhat - r));
                x = obj.thetaC(end-i+1, j);
            elseif z > obj.ztop
                % pull from thetaT
                [~, i] = min(abs(obj.zhat + obj.ztop - z));
                [~, j] = min(abs(obj.rtop - r));
                x = obj.thetaT(i, j);
            elseif z > obj.ztop + obj.h
                % set equal to ambient temp inside storage bin
                x = obj.thetaA;
            else
                % pull from thetaS
                [~, i] = min(abs(obj.zbar - z));
                [~, j] = min(abs(obj.rbar - r));
                x = obj.thetaS(i, j);
            end
        end
        function assembleUz(obj)
            % sets z-velocity component for every cell in mesh
            if isempty(obj.wbar)
                computeWbar(obj);
            end
            obj.uz = cell(7, 1);
            % stagnant region
            zs = 0:dz_:obj.ztop; rs = obj.a0:dr_:obj.b;
            obj.uz{2} = zeros(length(zs), length(rs));
            % top flowing surface
            zt = obj.ztop:dz_:obj.ztop+obj.h; rt = obj.a:dr_:obj.b;
            obj.uz{3} = zeros(length(zt)-1, length(rt));
            % center flowing channel
            zc = 0:dz_:obj.ztop; rc = 1e-6:dr_:obj.a0;
            wbar_ = interp1(obj.zcenter, obj.wbar, zc)';
            obj.uz{4} = repmat(wbar_, 1, length(rc)-1);
            % mixing region
            zm = obj.ztop:dz_:obj.ztop+obj.h; rm = 1e-6:dr_:obj.a;
            obj.uz{5} = zeros(length(zm)-1, length(rm)-1);
            % full domain
            obj.uz{1} = [obj.uz{5}, obj.uz{3}; obj.uz{4}, obj.uz{2}];
            obj.uz{6} = [zs, zt(2:end)]; obj.uz{7} = [rc(1:end-1), rs];
        end
        function assembleUr(obj)
            % sets r-velocity component for every cell in mesh
            if isempty(obj.ubar)
                computeUbar(obj);
            end
            obj.ur = cell(7, 1);
            dz_ = obj.dz; dr_ = obj.dr;
            % stagnant region
            zs = 0:dz_:obj.ztop; rs = obj.a0:dr_:obj.b;
            obj.ur{2} = zeros(length(zs), length(rs));
            % top flowing region
            zt = obj.ztop:dz_:obj.ztop+obj.h; rt = obj.a:dr_:obj.b;
            ubar_ = interp1(obj.rtop, obj.ubar, rt);
            obj.ur{3} = repmat(ubar_, length(zt)-1, 1);
            % center flowing channel
            zc = 0:dz_:obj.ztop; rc = 1e-6:dr_:obj.a0;
            obj.ur{4} = zeros(length(zc), length(rc)-1);
            % mixing region
            zm = obj.ztop:dz_:obj.ztop+obj.h; rm = 1e-6:dr_:obj.a;
            obj.ur{5} = zeros(length(zm)-1, length(rm)-1);
            % full domain
            obj.ur{1} = [obj.ur{5}, obj.ur{3}; obj.ur{4}, obj.ur{2}]; 
        end
        function assembleEnergyD(obj)
            % sets non-dimensional temperature values and velocities 
            % in a matching mesh for computing the energy equation in 
            % discharge mode
            if isempty(obj.thetaMatchD), obj.thetaMatchD = cell(8, 1); end
            if isempty(obj.ubar), computeUbar(obj); end
            if isempty(obj.wbar), computeWbar(obj); end
            obj.ur = cell(5, 1);
            obj.uz = cell(5, 1);
            dz_ = obj.dz; dr_ = obj.dr;
            % stagnant region
            zs = 0:dz_:obj.ztop; rs = obj.a0:dr_:obj.b;
            [Rs, Zs] = meshgrid(obj.rbar, obj.zbar);
            [Rsq, Zsq] = meshgrid(rs, zs);
            obj.thetaMatchD{5} = interp2(Rs, Zs, obj.thetaS, ...
                                                   Rsq, Zsq, 'makima');
            obj.ur{2} = zeros(length(zs), length(rs));
            obj.uz{2} = zeros(length(zs), length(rs));
            % top flowing region
            zt = obj.ztop:dz_:obj.ztop+obj.h; rt = obj.a:dr_:obj.b;
            [Rt, Zt] = meshgrid(obj.rtop, obj.zhat);
            [Rtq, Ztq] = meshgrid(rt, zt(2:end));
            obj.thetaMatchD{6} = interp2(Rt, Zt, obj.thetaT, ...
                                                   Rtq, Ztq, 'makima');
            ubar_ = interp1(obj.rtop, obj.ubar, rt);
            obj.ur{3} = repmat(ubar_, length(zt)-1, 1);
            obj.uz{3} = zeros(length(zt)-1, length(rt));
            % center flowing channel
            zc = 0:dz_:obj.ztop; rc = 1e-6:dr_:obj.a0;
            [Rc, Zc] = meshgrid(obj.rhat, obj.zcenter);
            [Rcq, Zcq] = meshgrid(rc(1:end-1), zc);
            obj.thetaMatchD{7} = interp2(Rc, Zc, obj.thetaC, ...
                                                 Rcq, Zcq, 'makima');
            obj.ur{4} = zeros(length(zc), length(rc)-1);
            wbar_ = interp1(obj.zcenter, obj.wbar, zc)';
            obj.uz{4} = repmat(wbar_, 1, length(rc)-1);
            % mixing region
            zm = obj.ztop:dz_:obj.ztop+obj.h; rm = 1e-6:dr_:obj.a;
            [Rmq, Zmq] = meshgrid(rm(1:end-1), zm(2:end));
            obj.thetaMatchD{8} = obj.thetaChat*ones(length(zm(1:end-1)), ...
                                                        length(rm(2:end)));
            obj.ur{5} = zeros(length(zm)-1, length(rm)-1);
            obj.uz{5} = zeros(length(zm)-1, length(rm)-1);
            % full domain
            Z = obj.thetaMatchD{3}; R = obj.thetaMatchD{4};
            Zq = [Zcq, Zsq; Zmq, Ztq]; Rq = [Rcq, Rsq; Rmq, Rtq];
            if isempty(obj.thetaMatchD{2})                
                obj.thetaMatchD{1} = [obj.thetaMatchD{7}, obj.thetaMatchD{5}; ...
                                     obj.thetaMatchD{8}, obj.thetaMatchD{6}]; 
                obj.thetaMatchD{2} = obj.thetaMatchD{1};
            else
                obj.thetaMatchD{2} = interp2(R, Z, obj.thetaMatchD{1}, ...
                                                        Rq, Zq, 'makima');
                obj.thetaMatchD{1} = [obj.thetaMatchD{7}, obj.thetaMatchD{5}; ...
                                     obj.thetaMatchD{8}, obj.thetaMatchD{6}];             
            end 
            obj.ur{1} = [obj.ur{4}, obj.ur{2}; obj.ur{5}, obj.ur{3}];
            obj.uz{1} = [obj.uz{4}, obj.uz{2}; obj.uz{5}, obj.uz{3}];
            obj.thetaMatchD{3} = Zq; obj.thetaMatchD{4} = Rq;
        end
        function assembleEnergyH(obj)
            % sets non-dimensional temperature values in a matching mesh 
            % for computing the energy equation in charge and holding mode
            if isempty(obj.thetaMatchH), obj.thetaMatchH = cell(10, 1); end
            zp = 0:obj.dz:obj.ztop; rp = 1e-6:obj.dr:obj.b;
            [Rp, Zp] = meshgrid(obj.rbarH, obj.zbarH);
            [Rpq, Zpq] = meshgrid(rp, zp);
            t = interp2(Rp, Zp, full(obj.thetaH), Rpq, Zpq);
            if isempty(obj.thetaMatchH{2})
                obj.thetaMatchH{1} = t; obj.thetaMatchH{2} = t;
            else
                Z = obj.thetaMatchH{3}; R = obj.thetaMatchH{4}; 
                obj.thetaMatchH{2} = interp2(R, Z, obj.thetaMatchH{1}, ...
                                                      Rpq, Zpq, 'makima');
                obj.thetaMatchH{1} = t;
            end
            obj.thetaMatchH{3} = Zpq; obj.thetaMatchH{4} = Rpq;  
        end
        function [theta_, z_, r_] = assembleTheta(obj)
            % sets non-dimensional temperature values in a matching mesh
            % set matching mesh size as smallest dz, dr
%             dz_ = min([obj.dzbar, obj.dzc, obj.dzhat]);
            dz_ = min([obj.dzc, obj.dzhat]);
%             dr_ = min([obj.drbar, obj.drhat, obj.drtop]);
            dr_ = min([obj.drhat, obj.drtop]);
            % stagnant region
            zs = 0:dz_:obj.ztop; rs = obj.a0:dr_:obj.b;
            [Rs, Zs] = meshgrid(obj.rbar, obj.zbar);
            [Rsq, Zsq] = meshgrid(rs, zs);
            thetaS_ = interp2(Rs, Zs, obj.thetaS, Rsq, Zsq);
            % top flowing region
            zt = obj.ztop:dz_:obj.ztop+obj.h; rt = obj.a:dr_:obj.b;
            [Rt, Zt] = meshgrid(obj.rtop, obj.zhat);
            [Rtq, Ztq] = meshgrid(rt, zt(2:end));
            thetaT_ = interp2(Rt, Zt, obj.thetaT, Rtq, Ztq);
            % center flowing channel
            zc = 0:dz_:obj.ztop; rc = 1e-6:dr_:obj.a0;
            [Rc, Zc] = meshgrid(obj.rhat, obj.zcenter);
            [Rcq, Zcq] = meshgrid(rc(1:end-1), zc);
            thetaC_ = interp2(Rc, Zc, obj.thetaC, Rcq, Zcq);
            % mixing region (set to zero to null cell-to-cell analysis
            zm = obj.ztop:dz_:obj.ztop+obj.h; rm = 1e-6:dr_:obj.a;
            thetaM_ = obj.thetaChat*ones(length(zm(1:end-1)), ...
                                                        length(rm(2:end)));
            % full domain
            theta_ = [thetaC_, thetaS_; thetaM_, thetaT_];
            z_ = [zs, zt]; r_ = [rc, rs];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % experimental validation and parameter fitting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function x = thetakExp(obj, k_)
            % returns theta matrix at time step k_
            ns = length(obj.zhat) - 1; nl = length(obj.zbar0); n = ns + nl;
            ms = length(obj.rhat) - 1; ml = length(obj.rbar); m = ms + ml;
            x = obj.thetaExp(n*(k_-1)+1:n*k_, m*(k_-1)+1:m*k_);
        end
        function addExperimentalData(obj, zbar_, rbar_, Fo_, theta_, convert)
            % note: length(Fo_) must equal length(theta_)
            % note: can only add for one location at a time
            ns = length(obj.zhat) - 1; nl = length(obj.zbar0); n = ns + nl;
            ms = length(obj.rhat) - 1; ml = length(obj.rbar); m = ms + ml;
            l = length(obj.thetaStoreExp);     
            obj.thetaStoreExp{l+1}.z = zbar_;
            obj.thetaStoreExp{l+1}.r = rbar_;
            obj.thetaStoreExp{l+1}.Fo = Fo_;
            obj.thetaStoreExp{l+1}.theta = theta_;
            if nargin > 4 && convert
                zbar_ = zbar_/obj.H;
                rbar_ = rbar_/obj.H;
                Fo_ = obj.t2Fo(Fo_);
                theta_ = obj.T2theta(theta_);
            end           
            % append experimental data set
            obj.zbarExp = [obj.zbarExp, zbar_];
            obj.rbarExp = [obj.rbarExp, rbar_];
            obj.FoExp = [obj.FoExp, Fo_];
            % form dimensionally matched experimental data matrix
            if isempty(obj.thetaExp)
                obj.thetaExp = sparse(kron(zeros(length(obj.Fo)), ...
                    zeros(n, m)));
            end
            % match experimental time, height and radial locations
            [~, zi] = min(abs(zbar_ - obj.z));
            obj.zbarExp_ = [obj.zbarExp_, obj.z(zi)];
            [~, ri] = min(abs(rbar_ - obj.r));
            obj.rbarExp_ = [obj.rbarExp_, obj.r(ri)];
            fi = zeros(size(Fo_));
            for i = 1:length(Fo_) 
                [~, fi(i)] = min(abs(double(obj.Fo) - double(Fo_(i))));
            end
            % place experimental temperature data points  
            for i = 1:length(fi)
                obj.thetaExp(length(obj.z)*(fi(i)-1)+zi, ...
                    length(obj.r)*(fi(i)-1)+ri) = theta_(i);
            end
            obj.thetaExp = sparse(obj.thetaExp);
        end
        function [result, error, runs] = fitWallConvection(obj, Z, R, ...
                                                               U, ShowPlot)
            % fits overall convection coefficients at side-wall and bottom
            % of storage bin
            % Z/R: z/r-locations to use for parameter fitting
            %      (Z and R need to be the same size and experimental data
            %      should be loaded close to these locations)
            % U: vector of overall heat transfer coefficients to sweep
            if nargin < 5
                ShowPlot = 1;
            end
            obj.zStore = Z; obj.rStore = R;
            result = struct('solution', [], 'thetaFit', [], ...
                            'thetaExp', [], ...
                            'rmse_total', [], 'U', []);
            runs{length(U)} = [];
            error{length(U)} = [];
            for i = 1:length(U)
                % simulate with sweeping h2 and h4 to minimize error
                obj.h2 = U(i); obj.h4 = U(i);                 
                resetPrimal(obj);
                reInitObj(obj);
                simulateDischargeTheta(obj, 'point');
                runs{i} = obj.thetaStore;
                error{i} = computeRMSE(obj);
                result.rmse_total = [result.rmse_total, ...
                        error{i}{1}.rmse];
                if length(error{i}) > 1 
                    for ii = 2:length(error{i})
                        result.rmse_total(i) = result.rmse_total(i) ...
                                              + error{i}{ii}.rmse;
                    end
                end        
            end
            [result.solution.rmse_total, idx] = min(result.rmse_total);
            result.solution.U = U(idx);
            result.thetaFit = runs{idx};
            result.thetaExp = obj.thetaStoreExp;
            result.U = U;
            if ShowPlot
                plotFitResults(obj, result);
            end                       
        end
        function [result, error, runs] = fitBoundaryConvection(obj, Z, ...
                                                            R, h, ShowPlot)
            % fits overall convection coefficients at side-wall and bottom
            % of storage bin
            % Z/R: z/r-locations to use for parameter fitting
            %      (Z and R need to be the same size and experimental data
            %      should be loaded close to these locations)
            % U: vector of heat transfer coefficients to sweep
            if nargin < 5
                ShowPlot = 1;
            end
            obj.zStore = Z; obj.rStore = R;
            result = struct('solution', [], 'thetaFit', [], ...
                            'thetaExp', [], ...
                            'rmse_total', [], 'h', []);
            runs{length(h)} = [];
            error{length(h)} = [];
            for i = 1:length(h)
                % simulate with sweeping h2 and h4 to minimize error 
                obj.a0 = h(i);
                resetPrimal(obj);
                reInitObj(obj);
                simulateDischargeTheta(obj, 'point');
                runs{i} = obj.thetaStore;
                error{i} = computeRMSE(obj);
                result.rmse_total = [result.rmse_total, ...
                        error{i}{1}.rmse];
                if length(error{i}) > 1 
                    for ii = 2:length(error{i})
                        result.rmse_total(i) = result.rmse_total(i) ...
                                              + error{i}{ii}.rmse;
                    end
                end        
            end
            [result.solution.rmse_total, idx] = min(result.rmse_total);
            result.solution.h = h(idx);
            result.thetaFit = runs{idx};
            result.thetaExp = obj.thetaStoreExp;
            result.h = h;
            if ShowPlot
                plotFitResults(obj, result);
            end                       
        end
        function e = computeRMSE(obj)
            % computes the rmse error between the computed thetaPoints and
            % the stored experimental data at matching points
            e{length(obj.zStore)} = [];
            for i = 1:length(e)
                ySim = interp1(obj.Fo, obj.thetaStore{i}.theta, ...
                    obj.FoExp);
                [~, ie] = min(abs(obj.zStore(i) - obj.zbarExp) ...
                        + abs(obj.rStore(i) - obj.rbarExp));
                yExp = obj.thetaStoreExp{ie}.theta;
                N = length(ySim);
                e{i} = struct('z', obj.zStore(i), 'r', obj.rStore(i), ...
                    'rmse', 1/sqrt(N)*norm(ySim - yExp));
            end            
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
        function y = S2T(obj)
            % matches size of a row of stagnant region with the the size
            % of a row in the top boundary
            n = ceil(size(obj.thetaS, 1)/2);
%             y = interp1(obj.rbar, obj.thetaS(n, :), obj.rtop, ...
%                  'linear', 'extrap');
            y = interp1(obj.rbar, obj.thetaS(end-1, :), obj.rtop, ...
                 'linear', 'extrap');
        end
        function y = S2C(obj)
            % matches size of a row of stagnant region with the the size
            % of a row in the top boundary
            y = interp1(obj.zbar, flipud(obj.thetaS(:, 1)), obj.zcenter, ...
                 'linear', 'extrap');
        end
        function y = T2S(obj, averaging)
            % matches size of a row in the top boundary with the size of a
            % row in the stagnant region
            if nargin < 1
                averaging = 'full';
            end
            if isempty(obj.thetaChat)
                computeThetaChat(obj);
            end
            switch averaging
                case 'full'
                    y = NaN*ones(size(obj.thetaT, 1), size(obj.thetaS, 2));
                    for i = 1:size(y, 1)
                        y(i, :) = interp1(obj.rtop, obj.thetaT(i, :), ...
                            obj.rbar, 'linear', 'extrap');
                    end
                case 'mid'
                    y = interp1(obj.rtop, mean(obj.thetaT), ...
                        obj.rbar, 'linear', 'extrap');
                case 'top'
                    y = interp1(obj.rtop, obj.thetaT(1, :), ...
                        obj.rbar, 'linear', 'extrap');
                case 'bottom'
                    y = interp1(obj.rtop, obj.thetaT(end, :), ...
                        obj.rbar, 'linear', 'extrap');
            end
        end
        function y = C2S(obj, averaging)
            % matches size of a column in the center boundary with the 
            % size of a column in the stagnant region
            if nargin < 1
                averaging = 'full';
            end
            switch averaging
                case 'full'
                    y = NaN*ones(size(obj.thetaS, 1), size(obj.thetaC, 2));
                    for j = 1:size(y, 2)
                            y(:, j) = interp1(obj.zcenter, ...
                                obj.thetaC(:, j), ...
                                obj.zbar, 'linear', 'extrap');
                    end
                case 'mid'
                    y = interp1(obj.zcenter, flipud(mean(obj.thetaC, 2)), ...
                        obj.zbar, 'linear', 'extrap');
                case 'top'
                    y = interp1(obj.zcenter, flipud(obj.thetaC(:, 1)), ...
                        obj.zbar, 'linear', 'extrap');
                case 'bottom'
                    y = interp1(obj.zcenter, flipud(obj.thetaC(:, end)), ...
                        obj.zbar, 'linear', 'extrap');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conservation validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function storeTotalMass(obj)
            % stores the actual mass contained in the tank at the current
            % time step
            if isempty(obj.simMass)
                obj.simMass = [];
            end
            obj.simMass = [obj.simMass; obj.FoNow, ...
                           obj.ztop*obj.H*pi*(obj.b*obj.H)^2*obj.rhoPack];
        end
        function computeTheoreticalMass(obj)
            % computes the theoretical mass contained in the tank over the
            % entire simulation period defined by Fo
            H_ = obj.H;
            obj.thryMass = [];
            for i = 1:length(obj.Fo)                
                obj.thryMass = [obj.thryMass; obj.Fo(i), ...
                                   H_*pi*(obj.b*obj.H)^2*obj.rhoPack];  
                H_ = H_ - obj.dH;
            end
        end
        function computeContinuity(obj)
            % computes the 2D continuity equation for every cell in the
            % domain. deviation from zero represents error (failing mass
            % conservation).
            [dz_, dr_] = assembleUz(obj); assembleUr(obj); 
            n = size(obj.ur{1}, 1); m = size(obj.ur{1}, 2);
            Dr = diffR(obj, n, m, dr_); Dz = diffZ(obj, n, m, dz_);
            mr1 = 1./obj.ur{7};
            mr2 = reshape(Dr*reshape((obj.ur{7}.*obj.ur{1})', [], 1), m, n)';
            mz = reshape(Dz*reshape(obj.uz{1}', [], 1), m, n)';
            obj.continuity = mr1 + mr2 + mz;                       
        end
        function computeEnergyD(obj)
            % computes the NS 2D energy equation at every cell in the
            % domain. deviation from zero represents error (failing energy
            % conservation).
            assembleUz(obj); assembleUr(obj); 
            assembleEnergyD(obj); t_ = obj.thetaMatchD;
            dz_ = obj.dz; dr_ = obj.dr;
            n = size(obj.ur{1}, 1); m = size(obj.ur{1}, 2);
            Dr = diffZ(obj, n, m, dr_); Dz = diffZ(obj, n, m, dz_);
            D2r = diff2R(obj, n, m, dr_); D2z = diff2Z(obj, n, m, dz_);
            ef = (t_{1} - t_{2})/obj.df;
            era = reshape(Dr*reshape((obj.ur{1}.*t_{2})', [], 1), m, n)';
            eza = reshape(Dz*reshape((obj.uz{1}.*t_{2})', [], 1), m, n)';
            erd1 = reshape(Dr*reshape((t_{2})', [], 1), m, n)'./t_{4};
            erd2 = reshape(D2r*reshape((t_{2})', [], 1), m, n)';
            ezd2 = reshape(D2z*reshape((t_{2})', [], 1), m, n)';
            obj.energy = ef + obj.Pe*(era + eza) - (erd1 + erd2 + ezd2);            
        end
        function computeEnergyH(obj)
            % computes heat equation for every cell for holding and
            % charging modes. deviation from zero represents error (failing
            % energy conservation)
            assembleEnergyH(obj); t_ = obj.thetaMatchH;
            dz_ = obj.dz; dr_ = obj.dr;
            n = size(t_{1}, 1); m = size(t_{1}, 2);
            Dr = diffR(obj, n, m, dr_);
            D2r = diff2R(obj, n, m, dr_); D2z = diff2Z(obj, n, m, dz_);
            ef = (t_{1} - t_{2})/obj.df;
            erd1 = reshape(Dr*reshape((t_{2})', [], 1), m, n)'./t_{4};
            erd2 = reshape(D2r*reshape((t_{2})', [], 1), m, n)';
            ezd2 = reshape(D2z*reshape((t_{2})', [], 1), m, n)';
            obj.energy = ef - (erd1 + erd2 + ezd2);
        end
        function u = computeInternalEnergy(obj, z_, r_, t)
            % computes the total internal energy contained in cylindrical
            % portion of storage bin for a given temperature distribution,
            % t with a temperature-dependent specific heat. t must be in
            % Kelvin
            c = obj.variableC(t);
            fr = simpsonIntegrator(obj, r_);
            fz = simpsonIntegrator(obj, z_);
            Ir = 2*pi*obj.rhoPack*(t.*r_.*c)*fr';
            u = Ir'*fz';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % accuracy testing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        function e = checkStagnantIntegralAccuracy(obj, convergence, ...
                                                          nrbar_, nzbar_)
            % tests the accuracy of the integral approximation used for
            % Fourier coefficients throughout the simulation by comparing
            % the same method to the analytical solution with a unity
            % temperature profile (approximate accuracy prediction)
            obj.ztop = 1;
            % set domain non-dimensional temperature to unity to allow for
            % analytical computation            
            % compute approximate (thetaIC) and analytical (Psi) solutions
            % and evaluate error with matrix norm
            if nargin < 2 || convergence == 0
                resetPrimal(obj);
                computeThetaIC(obj); thetaIC_ = obj.thetaIC; 
                resetPrimal(obj);
                computePsi(obj); 
                e = max(max(abs(thetaIC_ - obj.psiS)./obj.psiS))*100;
            else
                if nargin < 4
                    nzbar_ = logspace(1, 3, 10); 
                end
                if nargin < 3
                    nrbar_ = logspace(1, 3, 10);
                end
                e = NaN*ones(length(nrbar_), 1);
                for i = 1:length(nrbar_)
%                     obj.drbar = nrbar_(i);
%                     obj.dzbar = nzbar_(i);
                    obj.nrbar = ceil(nrbar_(i)); 
                    obj.nzbar = ceil(nzbar_(i));
                    reInitObj(obj); 
                    resetPrimal(obj); 
                    computeGCnm(obj, 0);
                    computeThetaIC(obj); thetaIC_ = obj.thetaIC; 
                    resetPrimal(obj);
                    computePsi(obj);
                    e(i) = max(max(abs(thetaIC_ - obj.psiS) ...
                                                        ./obj.psiS))*100;
                end
                figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.3], 'Visible', 'on');    
                hold on
                plot(nrbar_, e, '-k');
                title('Stagnant Integral Error', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylabel('error ($\%$)', 'interpreter', 'latex', ...
                    'FontSize', 14);
                xlabel('$dR$', 'interpreter', 'latex', 'FontSize', 14);
                set(gca, 'TickLabelInterpreter', 'latex')
                set(gcf, 'Color', [1 1 1])
%                 set(gca, 'XScale', 'log')
            end      
            % reinitialize primitive variables
            resetPrimal(obj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotting and vizualization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function s = printStatus(obj, theta_)
            % prints the current status of storage bin simulation
            if nargin < 2, theta_ = NaN; end
            s1 = sprintf('----------------------------------------- \n model time = %1.0f s \n', obj.Fo2t(obj.FoNow));
            s2 = sprintf(' prototype time = %1.0f s \n', obj.Fo2t(obj.FoNow, 1));
            s3 = sprintf(' storage mode = %s \n', obj.FoModeNow);
            if ~isnan(theta_)
                s4 = sprintf( ' max temperature = %1.0f °C \n', ...
                                              obj.theta2T(max(theta_, [], 'all')));
                s5 = sprintf( ' min temperature = %1.0f °C \n', ...
                                              obj.theta2T(min(theta_, [], 'all')));
            else
                s4 = '';
                s5 = '';
            end
            s6 = sprintf(' particle top = %1.2f \n----------------------------------------- \n', obj.ztop);  
            s = {s1; s2; s3; s4; s5; s6};
            display(s);
        end
        function animatePsi(obj, rr, filename)
            % compute theta if empty
            if isempty(obj.psiS)
                computePsi(obj, 1);
            end
            % initialize figure at Fo(1)
            n = length(obj.zbar); m = length(obj.rbar);
            figure('Units', 'normalized', ...
                'Position', [0 0 0.4*obj.b 0.4], 'Visible', 'off');            
            [~, pzrf] = plotZRPsi(obj, obj.psiS(1:n, 1:m), true);
            ylabel(colorbar, '$\Psi$', 'interpreter', 'latex', 'FontSize', 24);
            gif(filename, 'frame', gcf);
            % update sequentially
            for ii = 1:length(obj.Fo)
                gif;
                pzrf.ZData = obj.psiS((n*(ii-1)+1):n*ii, (m*(ii-1)+1):m*ii); 
                title(sprintf('$Fo$ = %1.2f', obj.Fo(ii)));
                pause(rr);
            end
        end
        function animateTheta(obj, rr, plot_filename, snapK, dimensions, ...
                                                                prototype)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ThetaPlot.gif'; end
            if nargin < 4, snapK = 1e10; end
            if nargin < 5, dimensions = 0; end
            if nargin < 6, prototype = 0; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = [obj.theta{2, 2}, ...
                              max(obj.theta{2, 2}), 1 + 1.000001*obj.h];
            r_ = obj.theta{2, 3};
            theta_  = [obj.theta{2, 1}; obj.thetaA*ones(2, length(r_))];
            zw_ = obj.zbarW; rw_ = []; thetaw_ = [];
            for i = 1:length(obj.rbarW)
                rw_ = [rw_, obj.rbarW{i}]; 
                thetaw_ = [thetaw_, obj.theta{2, 7}{i}];
            end
            if dimensions
                [~, tfigs1, tfigs2] = plotZRTemp(obj, obj.theta2T(theta_), ...
                                          z_, r_, obj.theta2T(thetaw_), ...
                                                          zw_, rw_, true);
                ylabel(colorbar, '$T$ ($^\circ$C)', 'interpreter', ...
                    'latex', 'FontSize', 14); 
                caxis([100, 800]);
            else
                [~, tfigs1, tfigs2] = plotZRTemp(obj, theta_, z_, r_, thetaw_, ...
                                                        zw_, rw_, true);
                ylabel(colorbar, '$\theta$', 'interpreter', 'latex', ...
                'FontSize', 14); 
                caxis([obj.thetaA, 1]);
            end
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 2:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = [obj.theta{ii_, 2}, ...
                              max(obj.theta{ii_, 2}), 1 + 1.000001*obj.h];
                r_ = obj.theta{ii_, 3};
                theta_  = [obj.theta{ii_, 1}; ...
                                  obj.thetaA*ones(2, length(r_))];                              
                [R, Z] = meshgrid(r_, z_); 
                thetaw_ = [];
                for i = 1:length(obj.rbarW)
                    thetaw_ = [thetaw_, obj.theta{ii_, 7}{i}];
                end
                [RW, ZW] = meshgrid(rw_, zw_);                
                if dimensions 
                    set(tfigs1, 'XData', R, 'YData', Z, ...
                      'ZData', obj.theta2T(theta_));
                    set(tfigs2, 'XData', RW, 'YData', ZW, ...
                      'ZData', obj.theta2T(thetaw_));
                else
                    set(tfigs1, 'XData', R, 'YData', Z, ...
                      'ZData', theta_);
                    set(tfigs2, 'XData', RW, 'YData', ZW, ...
                      'ZData', thetaw_);
                end
                title(sprintf('$t$ = %1.2f h', t/3600), 'interpreter', 'latex', ...
                            'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureDimensionalTemp(obj, ii_);
                end               
            end
        end
        function animateThetaNoWall(obj, rr, plot_filename, snapK, dimensions, ...
                                                                prototype)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ThetaPlot.gif'; end
            if nargin < 4, snapK = 1e10; end
            if nargin < 5, dimensions = 0; end
            if nargin < 6, prototype = 0; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = [obj.theta{2, 2}, ...
                              max(obj.theta{2, 2}), 1 + 1.000001*obj.h];
            r_ = [-fliplr(obj.theta{2, 3}), obj.theta{2, 3}];
            thetaA_ = obj.theta{2, 12};
            theta_  = [fliplr(obj.theta{2, 1}), obj.theta{2, 1}; ...
                           thetaA_*ones(2, length(r_))];
            if dimensions
                [~, tfigs1] = plotZRTempNoWall(obj, obj.theta2T(theta_), ...
                                          z_, r_, true);
                ylabel(colorbar, '$T$ ($^\circ$C)', 'interpreter', ...
                    'latex', 'FontSize', 14); 
                caxis([20, 800]);
            else
                [~, tfigs1] = plotZRTempNoWall(obj, theta_, z_, r_, true);
                ylabel(colorbar, '$\theta$', 'interpreter', 'latex', ...
                'FontSize', 14); 
                caxis([obj.thetaA, 1]);
            end
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 2:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = [obj.theta{ii_, 2}, ...
                              max(obj.theta{ii_, 2}), 1 + 1.000001*obj.h];
                r_ = [-fliplr(obj.theta{ii_, 3}), obj.theta{ii_, 3}];
                thetaA_ = obj.theta{ii_, 12};
                theta_  = [fliplr(obj.theta{ii_, 1}), obj.theta{ii_, 1}; ...
                                  thetaA_*ones(2, length(r_))];                              
                [R, Z] = meshgrid(r_, z_);                 
                if dimensions 
                    set(tfigs1, 'XData', R, 'YData', Z, ...
                      'ZData', obj.theta2T(theta_));
                else
                    set(tfigs1, 'XData', R, 'YData', Z, ...
                      'ZData', theta_);
                end
                title(sprintf('$t$ = %1.2f h', t/3600), 'interpreter', 'latex', ...
                            'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureDimensionalTemp(obj, ii_);
                end               
            end
        end
        function animateThetaK(obj, rr, plot_filename, snapK, dimensions, ...
                                                                prototype)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ThetaPlot.gif'; end
            if nargin < 4, snapK = 1e10; end
            if nargin < 5, dimensions = 0; end
            if nargin < 6, prototype = 0; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadThetaK(obj, length(obj.Fo));  
            else, loadThetaK(obj, obj.ls); end
            z_ = [obj.thetaK{2, 2}, ...
                              max(obj.thetaK{2, 2}), 1 + 1.000001*obj.h];
            r_ = obj.thetaK{2, 3};
            theta_  = [obj.thetaK{2, 1}; obj.thetaA*ones(2, length(r_))];
            zw_ = obj.zbarW; rw_ = []; thetaw_ = [];
            for i = 1:length(obj.rbarW)
                rw_ = [rw_, obj.rbarW{i}]; 
                thetaw_ = [thetaw_, obj.thetaK{2, 6}{i}];
            end
            if dimensions
                [~, tfigs1, tfigs2] = plotZRTemp(obj, obj.theta2T(theta_), ...
                                          z_, r_, obj.theta2T(thetaw_), ...
                                                          zw_, rw_, true);
                ylabel(colorbar, '$T$ ($^\circ$C)', 'interpreter', ...
                    'latex', 'FontSize', 14); 
                caxis([100, 800]);
            else
                [~, tfigs1, tfigs2] = plotZRTemp(obj, theta_, z_, r_, thetaw_, ...
                                                        zw_, rw_, true);
                ylabel(colorbar, '$\theta$', 'interpreter', 'latex', ...
                'FontSize', 14); 
                caxis([obj.thetaA, 1]);
            end
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 2:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadThetaK(obj, ii+obj.ls-1);
                end
                z_ = [obj.thetaK{ii_, 2}, ...
                              max(obj.thetaK{ii_, 2}), 1 + 1.000001*obj.h];
                r_ = obj.thetaK{ii_, 3};
                theta_  = [obj.thetaK{ii_, 1}; ...
                                  obj.thetaA*ones(2, length(r_))];                              
                [R, Z] = meshgrid(r_, z_); 
                thetaw_ = [];
                for i = 1:length(obj.rbarW)
                    thetaw_ = [thetaw_, obj.thetaK{ii_, 6}{i}];
                end
                [RW, ZW] = meshgrid(rw_, zw_);                
                if dimensions 
                    set(tfigs1, 'XData', R, 'YData', Z, ...
                      'ZData', obj.theta2T(theta_));
                    set(tfigs2, 'XData', RW, 'YData', ZW, ...
                      'ZData', obj.theta2T(thetaw_));
                else
                    set(tfigs1, 'XData', R, 'YData', Z, ...
                      'ZData', theta_);
                    set(tfigs2, 'XData', RW, 'YData', ZW, ...
                      'ZData', thetaw_);
                end
                title(sprintf('$t$ = %1.2f s', t), 'interpreter', 'latex', ...
                            'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureDimensionalTemp(obj, ii_);
                end               
            end
        end
        function animateWallLineK(obj, rr, plot_filename, tEnd)            
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'WallLine.gif'; end  
            if nargin < 4, tEnd = obj.Fo2t(obj.Fo(end), 1); end
            % wall temp animation for t1
            figure('Units', 'normalized', 'color', 'white', ...
                'Position', [0 0 0.5 0.3]);  
            if length(obj.Fo) < obj.ls, loadThetaK(obj, length(obj.Fo));  
            else, loadThetaK(obj, obj.ls); end
            zw_ = objzbar.W; rw_ = []; thetaw_ = [];
            [~, zwi_] = min(abs(zw_ - 0.5));
            for i = 1:length(obj.rbarW)
                rw_ = [rw_, obj.rbarW{i}*obj.Hp]; 
                thetaw_ = [thetaw_, obj.thetaK{1, 6}{i}(zwi_, :)];
            end
            l1 = plot(rw_, obj.theta2T(thetaw_), '-k', 'LineWidth', 2);
            xlim([min(rw_), max(rw_)])
            ylim([0, 100])
            timeTitle = title( ...
                   sprintf('$t$ = %1.0f h', 0), 'Interpreter', 'latex', 'FontSize', 14);
            xlabel('$r$ ($m$)', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$T$ ($^\circ C$)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 12)
            [~, k_end] = min(abs(tEnd - obj.Fo2t(obj.Fo, 1)));
            gif(plot_filename, 'frame', gcf);
            for k_ = 2:k_end
                gif;
                kk_ = mod(k_-1, obj.ls) + 1;
                t_ = obj.Fo2t(obj.Fo(k_), 1);
                if kk_ == 1 && k_ ~= length(obj.Fo)
                    loadThetaK(obj, k_+obj.ls-1);
                end
                thetaw_ = [];
                for i = 1:length(obj.rbarW)
                    thetaw_ = [thetaw_, obj.thetaK{kk_, 6}{i}(zwi_, :)];
                end
                timeTitle.String = sprintf('$t$ = %1.0f h', t_);
                set(l1, 'YData', obj.theta2T(thetaw_));
                pause(rr);
            end
        end
        function animateWallLine(obj, rr, plot_filename, tEnd)            
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'WallLine.gif'; end  
            if nargin < 4, tEnd = obj.Fo2t(obj.Fo(end), 1); end
            % wall temp animation for t1
            figure('Units', 'normalized', 'color', 'white', ...
                'Position', [0 0 0.5 0.3]);  
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            zw_ = obj.zbarW; rw_ = []; thetaw_ = [];
            [~, zwi_] = min(abs(zw_ - 0.5));
            for i = 1:length(obj.rbarW)
                rw_ = [rw_, obj.rbarW{i}*obj.Hp]; 
                thetaw_ = [thetaw_, obj.theta{1, 7}{i}(zwi_, :)];
            end
            l1 = plot(rw_, obj.theta2T(thetaw_), '-k', 'LineWidth', 2);
            xlim([min(rw_), max(rw_)])
            ylim([0, 100])
            timeTitle = title( ...
                   sprintf('$t$ = %1.0f h', 0), 'Interpreter', 'latex', 'FontSize', 14);
            xlabel('$r$ ($m$)', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$T$ ($^\circ C$)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...l
                'TickLabelInterpreter', 'latex', 'FontSize', 12)
            [~, k_end] = min(abs(tEnd - obj.Fo2t(obj.Fo, 1)));
            gif(plot_filename, 'frame', gcf);
            for k_ = 2:k_end
                gif;
                kk_ = mod(k_-1, obj.ls) + 1;
                t_ = obj.Fo2t(obj.Fo(k_), 1);
                if kk_ == 1 && k_ ~= length(obj.Fo)
                    loadTheta(obj, k_+obj.ls-1);
                end
                thetaw_ = [];
                for i = 1:length(obj.rbarW)
                    thetaw_ = [thetaw_, obj.theta{kk_, 6}{i}(zwi_, :)];
                end
                timeTitle.String = sprintf('$t$ = %1.0f h', t_/3600);
                set(l1, 'YData', obj.theta2T(thetaw_));
                pause(rr);
            end
        end
        function animateCenterTheta(obj, rr, plot_filename)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'CenterTheta.gif'; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = obj.theta{2, 2};
            [~, i] = min(abs((z_(end) - obj.h) - z_)); z_ = z_(1:i);           
            r_ = obj.rhat;
            theta_ = obj.theta2T(obj.theta{2, 1}(1:i, 1:length(obj.rhat)));
            [~, thetafigs] = plotZRCenterSensitivity(obj, theta_, z_, r_, true);
            ylabel(colorbar, '$T$ ($^\circ$C)', ...
                'interpreter', 'latex', 'FontSize', 14); 
            colormap(jet);
            caxis([600, max(obj.T0(:))]);
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 3:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii));
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = obj.theta{ii_, 2};
                [~, i] = min(abs((z_(end) - obj.h) - z_)); z_ = z_(1:i);              
                r_ = obj.rhat;
                theta_ = obj.theta2T(obj.theta{ii_, 1}(1:i, 1:length(obj.rhat)));
                [R_, Z_] = meshgrid(r_, z_); 
                [R, Z] = meshgrid(r_, linspace(0, z_(end), size(theta_, 1)));
                theta_ = interp2(R, Z, theta_, R_, Z_);
                Z_ = [Z_; ones(1, length(r_))*max(z_); ...
                        ones(1, length(r_)) + 1e-6];
                R_ = [R_; R_(1, :); R_(1, :)];
                theta_  = [theta_; obj.thetaA*ones(2, length(r_))];
                set(thetafigs, 'XData', R_, 'YData', Z_, 'ZData', ... 
                    theta_);
                ttext = '$t$ = %1.2f s';
                title(sprintf(ttext, t), ...
                               'interpreter', 'latex', 'FontSize', 14);
                pause(rr);              
            end
        end
        function animateCenterScPe(obj, rr, plot_filename)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ScPeGif.gif'; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = obj.theta{2, 2};
            [~, i] = min(abs((z_(end) - obj.h) - z_)); z_ = z_(1:i);           
            r_ = obj.rhat;
            scPe_  = flipud(obj.theta{2, 10});
            [~, scPefigs] = plotZRCenterSensitivity(obj, scPe_, z_, r_, true);
            ylabel(colorbar, '$\frac{1}{\overline{s}_{Pe}}\frac{\partial \theta}{\partial Pe}$', ...
                'interpreter', 'latex', 'FontSize', 14); 
            caxis([0, 1]);
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 3:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii));
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = obj.theta{ii_, 2};
                [~, i] = min(abs((z_(end) - obj.h) - z_)); z_ = z_(1:i);              
                r_ = obj.rhat;
                scPe_  = flipud(obj.theta{ii_, 10});
                [R_, Z_] = meshgrid(r_, z_); 
                [R, Z] = meshgrid(r_, linspace(0, z_(end), size(scPe_, 1)));
                scPe_ = interp2(R, Z, scPe_, R_, Z_);
                Z_ = [Z_; ones(1, length(r_))*max(z_); ...
                        ones(1, length(r_)) + 1e-6];
                R_ = [R_; R_(1, :); R_(1, :)];
                scPe_  = [scPe_; max(scPe_(:))*ones(2, length(r_))];
                set(scPefigs, 'XData', R_, 'YData', Z_, 'ZData', ... 
                    scPe_/max(scPe_(:)));
                ttext = '$t$ = %1.2f s \n $\\overline{s}_{Pe}$ = %1.2e';
                title(sprintf(ttext, t, max(scPe_(:))), ...
                               'interpreter', 'latex', 'FontSize', 14);
                pause(rr);              
            end
        end
        function animateCenterScQc(obj, rr, plot_filename)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ScQcGif.gif'; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = obj.theta{2, 2};
            [~, i] = min(abs((z_(end) - obj.h) - z_)); z_ = z_(1:i);           
            r_ = obj.rhat;
            scQc_  = flipud(obj.theta{2, 11});
            [~, scQcfigs] = plotZRCenterSensitivity(obj, scQc_, z_, r_, true);
            ylabel(colorbar, '$\frac{1}{\overline{s}_{Q}}\frac{\partial \theta}{\partial Q}$', ...
                'interpreter', 'latex', 'FontSize', 14); 
            caxis([0, 1]);
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 3:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii));
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = obj.theta{ii_, 2};
                [~, i] = min(abs((z_(end) - obj.h) - z_)); z_ = z_(1:i);              
                r_ = obj.rhat;
                scQc_  = flipud(obj.theta{ii_, 11});
                [R_, Z_] = meshgrid(r_, z_); 
                [R, Z] = meshgrid(r_, linspace(0, z_(end), size(scQc_, 1)));
                scQc_ = interp2(R, Z, scQc_, R_, Z_);
                Z_ = [Z_; ones(1, length(r_))*max(z_); ...
                        ones(1, length(r_)) + 1e-6];
                R_ = [R_; R_(1, :); R_(1, :)];
                scQc_  = [scQc_; max(scQc_(:))*ones(2, length(r_))];
                set(scQcfigs, 'XData', R_, 'YData', Z_, 'ZData', ... 
                    scQc_/max(scQc_(:)));
                ttext = '$t$ = %1.2f s \n $\\overline{s}_{Q}$ = %1.2e';
                title(sprintf(ttext, t, max(scQc_(:))), ...
                               'interpreter', 'latex', 'FontSize', 14);
                pause(rr);              
            end
        end
        function animateCenterSca(obj, rr, plot_filename)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ScaGif.gif'; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = obj.theta{2, 2};
            [~, i] = min(abs((z_(end) - obj.h) - z_)); z_ = z_(1:i);           
            r_ = obj.rhat;
            sca_  = flipud(obj.theta{2, 12});
            [~, scafigs] = plotZRCenterSensitivity(obj, sca_, z_, r_, true);
            ylabel(colorbar, '$\frac{1}{\overline{s}_{a}}\frac{\partial \theta}{\partial a}$', ...
                'interpreter', 'latex', 'FontSize', 14); 
            caxis([0, 1]);
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 3:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii));
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = obj.theta{ii_, 2};
                [~, i] = min(abs((z_(end) - obj.h) - z_)); z_ = z_(1:i);              
                r_ = obj.rhat;
                sca_  = flipud(obj.theta{ii_, 12});
                [R_, Z_] = meshgrid(r_, z_); 
                [R, Z] = meshgrid(r_, linspace(0, z_(end), size(sca_, 1)));
                sca_ = interp2(R, Z, sca_, R_, Z_);
                Z_ = [Z_; ones(1, length(r_))*max(z_); ...
                        ones(1, length(r_)) + 1e-6];
                R_ = [R_; R_(1, :); R_(1, :)];
                sca_  = [sca_; max(sca_(:))*ones(2, length(r_))];
                set(scafigs, 'XData', R_, 'YData', Z_, 'ZData', ... 
                    sca_/max(sca_(:)));
                ttext = '$t$ = %1.2f s \n $\\overline{s}_{a}$ = %1.2e';
                title(sprintf(ttext, t, max(sca_(:))), ...
                               'interpreter', 'latex', 'FontSize', 14);
                pause(rr);              
            end
        end
        function animateContinuity(obj, rr, plot_filename, snapK)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ContinuityPlot.gif'; end
            if nargin < 4, snapK = 1e10; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = [obj.theta{1, 8}, ...
                              max(obj.theta{1, 8}), 1 + 1.000001*obj.h];
            r_ = obj.theta{1, 9};
            c_  = [flipud(obj.theta{1, 6}); ones(2, length(r_))];
            [~, tfigs] = plotZRConservation(obj, abs(c_) + 1e-6, z_, r_, true);
            ylabel(colorbar, ...
                '$\frac{1}{r}\frac{\partial}{\partial r}\left(r\:u_r\right) + \frac{\partial u_z}{\partial z}$', ...
                'interpreter', 'latex', 'FontSize', 14); 
            caxis([1e-6, 1]);
            set(gca, 'ColorScale', 'log');
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 2:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii));
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = [obj.theta{ii_, 8}, ...
                              max(obj.theta{ii_, 8}), 1 + 1.000001*obj.h];
                r_ = obj.theta{ii_, 9};
                c_  = [flipud(obj.theta{ii_, 6}); ones(2, length(r_))];
                [R, Z] = meshgrid(r_, z_); 
                set(tfigs, 'XData', R, 'YData', Z, 'ZData', abs(c_) + 1e-6);
                title(sprintf('$t$ = %1.2f s', t), 'interpreter', 'latex', ...
                            'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureDimensionalTemp(obj, ii_);
                end               
            end
        end
        function animateEnergy(obj, rr, plot_filename, snapK)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'EnergyPlot.gif'; end
            if nargin < 4, snapK = 1e10; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = [obj.theta{1, 7}, ...
                              max(obj.theta{1, 7}), 1 + 1.000001*obj.h];
            r_ = obj.theta{1, 8};
            c_  = [flipud(obj.theta{1, 6}); ones(2, length(r_))];
            [~, tfigs] = plotZRConservation(obj, abs(c_) + 1e-6, z_, r_, true);
            ylabel(colorbar, ...
                '$\frac{\partial \theta}{\partial t} + Pe\left(\frac{\partial(u_r\:\theta)}{\partial r} + \frac{\partial(u_z\:\theta)}{\partial z}\right) - \left(\frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial \theta}{\partial r}\right) + \frac{\partial^2\theta}{\partial z^2}\right)$', ...
                'interpreter', 'latex', 'FontSize', 14); 
            caxis([1e-6, 1]);
%             set(gca, 'ColorScale', 'log');
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 2:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii));
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = [obj.theta{ii_, 7}, ...
                              max(obj.theta{ii_, 7}), 1 + 1.000001*obj.h];
                r_ = obj.theta{ii_, 8};
                c_  = [flipud(obj.theta{ii_, 6}); ones(2, length(r_))];
                [R, Z] = meshgrid(r_, z_); 
                set(tfigs, 'XData', R, 'YData', Z, 'ZData', abs(c_) + 1e-6);
                title(sprintf('$t$ = %1.2f s', t), 'interpreter', 'latex', ...
                            'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureDimensionalTemp(obj, ii_);
                end               
            end
        end
        function animateThetaH(obj, rr, plot_filename, snapK, dimensions)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ThetaH.gif'; end
            if nargin < 4, snapK = 1e10; end
            if nargin < 5, dimensions = 0; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = [obj.theta{1, 2}, ...
                              max(obj.theta{1, 2}), 1 + 1.000001*obj.h];
            r_ = obj.theta{1, 3};
            theta_  = [obj.theta{1, 1}; obj.thetaA*ones(2, length(r_))];
            if dimensions
                [~, tfigs] = plotZRTemp(obj, obj.theta2T(theta_), ...
                                                          z_, r_, true);
                ylabel(colorbar, '$T$ ($^\circ$C)', 'interpreter', ...
                    'latex', 'FontSize', 14); 
                caxis([500, max(obj.T0(:))]);
            else
                [~, tfigs] = plotZRTemp(obj, theta_, z_, r_, true);
                ylabel(colorbar, '$\theta$', 'interpreter', 'latex', ...
                'FontSize', 14); 
                caxis([obj.thetaA, 1]);
            end
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 2:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii));
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = [obj.theta{ii_, 2}, ...
                              max(obj.theta{ii_, 2}), 1 + 1.000001*obj.h];
                r_ = obj.theta{ii_, 3};
                theta_  = [obj.theta{ii_, 1}; ...
                                  obj.thetaA*ones(2, length(r_))];
                [R, Z] = meshgrid(r_, z_); 
                if dimensions 
                    set(tfigs, 'XData', R, 'YData', Z, ...
                        'ZData', obj.theta2T(theta_));
                else
                    tfigs.ZData = theta_; 
                    set(tfigs, 'XData', R, 'YData', Z, 'ZData', theta_);
                end
                title(sprintf('$t$ = %1.2f s', t), 'interpreter', 'latex', ...
                            'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureDimensionalTemp(obj, ii_);
                end               
            end
        end
        function animateThetaCh(obj, rr, plot_filename, snapK, dimensions)
            % generate full domain temperature gif to be saved
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'ThetaCh.gif'; end
            if nargin < 4, snapK = 1e10; end
            if nargin < 5, dimensions = 0; end
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = [obj.theta{1, 2}, ...
                              max(obj.theta{1, 2}), 1 + 1.000001*obj.h];
            r_ = obj.theta{1, 3};
            theta_  = [obj.theta{1, 1}; obj.thetaA*ones(2, length(r_))];
            if dimensions
                [~, tfigs] = plotZRTemp(obj, obj.theta2T(theta_), ...
                                                          z_, r_, true);
                ylabel(colorbar, '$T$ ($^\circ$C)', 'interpreter', ...
                    'latex', 'FontSize', 14); 
                caxis([500, max(obj.T0(:))]);
            else
                [~, tfigs] = plotZRTemp(obj, theta_, z_, r_, true);
                ylabel(colorbar, '$\theta$', 'interpreter', 'latex', ...
                'FontSize', 14); 
                caxis([obj.thetaA, 1]);
            end
            gif(plot_filename, 'frame', gcf);
            % update sequentially
            for ii = 2:length(obj.Fo)                                      
                gif;
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii));
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = [obj.theta{ii_, 2}, ...
                              max(obj.theta{ii_, 2}), 1 + 1.000001*obj.h];
                r_ = obj.theta{ii_, 3};
                theta_  = [obj.theta{ii_, 1}; ...
                                  obj.thetaA*ones(2, length(r_))];
                [R, Z] = meshgrid(r_, z_); 
                if dimensions 
                    set(tfigs, 'XData', R, 'YData', Z, ...
                        'ZData', obj.theta2T(theta_));
                else
                    tfigs.ZData = theta_; 
                    set(tfigs, 'XData', R, 'YData', Z, 'ZData', theta_);
                end
                title(sprintf('$t$ = %1.2f s', t), 'interpreter', 'latex', ...
                            'FontSize', 14);
                pause(rr);
                % save plot if on save time-step
                if mod(ii, snapK) == 0
                    captureDimensionalTemp(obj, ii_);
                end               
            end
        end
        function data = animateRadialTempLine(obj, rr, plot_filename, zCoord, dimensions)            
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'RadialLine.gif'; end  
            if nargin < 4, zCoord = 0.5; end
            if nargin < 5, dimensions = 1; end
            data = {};
            % wall temp animation for t1
            figure('Units', 'normalized', 'color', 'white', ...
                'Position', [0 0 0.5 0.3]);  
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = obj.theta{2, 2};
            r_ = obj.theta{2, 3};
            [~, zwi_] = min(abs(z_ - zCoord));
            if dimensions
                thetaLine = obj.theta2T(obj.theta{2, 1}(zwi_, :));
            else
                thetaLine = obj.theta{2, 1}(zwi_, :);
            end            
            l1 = plot(r_, thetaLine, '-k', 'LineWidth', 2);
            xlim([min(r_), max(r_)])           
            timeTitle = title( ...
                   sprintf('$t$ = %1.1f h', 0), 'Interpreter', 'latex', 'FontSize', 14);
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);  
            if dimensions
%                 ylim([600, 800])
                ylabel(sprintf('$T$($z/H$ = %1.2f) (degC)', zCoord), 'interpreter', 'latex', 'FontSize', 14);
            else
%                 ylim([0.6, 1])
                ylabel(sprintf('$\theta$($z/H$ = %1.2f)', zCoord), 'interpreter', 'latex', 'FontSize', 14);
            end
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 12)
            data{1, 1} = 0;     % (s) time
            data{1, 2} = r_;
            data{1, 3} = thetaLine;
            gif(plot_filename, 'frame', gcf);
            for k_ = 2:length(obj.Fo)
                gif;
                kk_ = mod(k_-1, obj.ls) + 1;
                t_ = obj.Fo2t(obj.Fo(k_), 1);
                if kk_ == 1 && k_ ~= length(obj.Fo)
                    loadTheta(obj, k_+obj.ls-1);
                end
                z_ = obj.theta{kk_, 2};
                r_ = obj.theta{kk_, 3};
                [~, zwi_] = min(abs(z_ - zCoord));
                if dimensions
                    thetaLine = obj.theta2T(obj.theta{kk_, 1}(zwi_, :));   
                else
                    thetaLine = obj.theta{kk_, 1}(zwi_, :);
                end
                timeTitle.String = sprintf('$t$ = %1.1f h', t_/3600);
                set(l1, 'XData', r_, 'YData', thetaLine);
                data{k_, 1} = t_;
                data{k_, 2} = r_;
                data{k_, 3} = thetaLine;
                pause(rr);
            end
        end
        function data = animateDiscreteParticleWallLine(obj, rr, plot_filename, zCoord, dimensions)
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'RadialLine.gif'; end  
            if nargin < 4, zCoord = 0.1; end
            if nargin < 5, dimensions = 1; end
            data = {};
            % wall temp animation for t1
            figure('Units', 'normalized', 'color', 'white', ...
                'Position', [0 0 0.5 0.3]);  
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = obj.theta{2, 2};
            r_ = [obj.theta{2, 3}*obj.Hp, obj.rWL(2:end)'];
            [~, zwi_] = min(abs(z_ - zCoord));
            if dimensions
                thetaLine = [obj.theta2T(obj.theta{2, 1}(zwi_, :)), ...
                    obj.theta2T(obj.theta{2, 7})];
            else
                thetaLine = [obj.theta{2, 1}(zwi_, :), obj.theta{2, 7}];
            end            
            l1 = plot(r_, thetaLine, '-k', 'LineWidth', 2);
            xlim([min(r_), max(r_)])           
            timeTitle = title( ...
                   sprintf('$t$ = %1.1f h', 0), 'Interpreter', 'latex', 'FontSize', 14);
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);  
            if dimensions
%                 ylim([600, 800])
                ylabel(sprintf('$T$($z/H$ = %1.2f) (degC)', zCoord), 'interpreter', 'latex', 'FontSize', 14);
            else
%                 ylim([0.6, 1])
                ylabel(sprintf('$\theta$($z/H$ = %1.2f)', zCoord), 'interpreter', 'latex', 'FontSize', 14);
            end
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 12)
            data{1, 1} = 0;     % (s) time
            data{1, 2} = r_;
            data{1, 3} = thetaLine;
            gif(plot_filename, 'frame', gcf);
            for k_ = 2:length(obj.Fo)
                gif;
                kk_ = mod(k_-1, obj.ls) + 1;
                t_ = obj.Fo2t(obj.Fo(k_), 1);
                if kk_ == 1 && k_ ~= length(obj.Fo)
                    loadTheta(obj, k_+obj.ls-1);
                end
                z_ = obj.theta{kk_, 2};
                r_ = [obj.theta{kk_, 3}*obj.Hp, obj.rWL(2:end)'];
                [~, zwi_] = min(abs(z_ - zCoord));
                if dimensions
                    thetaLine = [obj.theta2T(obj.theta{kk_, 1}(zwi_, :)), ...
                    obj.theta2T(obj.theta{kk_, 7})];   
                else
                    thetaLine = [obj.theta{kk_, 1}(zwi_, :), obj.theta{kk_, 7}];
                end
                timeTitle.String = sprintf('$t$ = %1.1f h', t_/3600);
                set(l1, 'XData', r_, 'YData', thetaLine);
                data{k_, 1} = t_;
                data{k_, 2} = r_;
                data{k_, 3} = thetaLine;
                pause(rr);
            end
        end  
        function data = animateVerticalTempLine(obj, rr, plot_filename, rCoord, dimensions)            
            if nargin < 2, rr = 0.5; end
            if nargin < 3, plot_filename = 'VerticalLine.gif'; end  
            if nargin < 4, rCoord = 0.1; end
            if nargin < 5, dimensions = 1; end
            data = {};
            % wall temp animation for t1
            figure('Units', 'normalized', 'color', 'white', ...
                'Position', [0 0 0.3 0.5]);  
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            z_ = obj.theta{2, 2};
            r_ = obj.theta{2, 3};
            [~, rwi_] = min(abs(r_ - rCoord));
            if dimensions
                thetaLine = obj.theta2T(obj.theta{2, 1}(:, rwi_));
            else
                thetaLine = obj.theta{2, 1}(:, rwi_);
            end            
            l1 = plot(thetaLine, z_, '-k', 'LineWidth', 2);
            ylim([min(z_), max(z_)])           
            timeTitle = title( ...
                   sprintf('$t$ = %1.1f h', 0), 'Interpreter', 'latex', 'FontSize', 14);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14); 
            if dimensions
%                 xlim([600, 800])
                xlabel(sprintf('$T$($r/H$ = %1.2f) (degC)', rCoord), 'interpreter', 'latex', 'FontSize', 14);
            else
%                 xlim([0.6, 1])
                xlabel(sprintf('$\theta$($r/H$ = %1.2f)', rCoord), 'interpreter', 'latex', 'FontSize', 14);
            end
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...l
                'TickLabelInterpreter', 'latex', 'FontSize', 12)
            data{1, 1} = 0;     % (s) time
            data{1, 2} = z_;
            data{1, 3} = thetaLine;
            gif(plot_filename, 'frame', gcf);
            for k_ = 2:length(obj.Fo)
                gif;
                kk_ = mod(k_-1, obj.ls) + 1;
                t_ = obj.Fo2t(obj.Fo(k_), 1);
                if kk_ == 1 && k_ ~= length(obj.Fo)
                    loadTheta(obj, k_+obj.ls-1);
                end
                z_ = obj.theta{kk_, 2};
                r_ = obj.theta{kk_, 3};
                [~, rwi_] = min(abs(r_ - rCoord));
                if dimensions
                    thetaLine = obj.theta2T(obj.theta{kk_, 1}(:, rwi_));   
                else
                    thetaLine = obj.theta{kk_, 1}(:, rwi_);
                end
                timeTitle.String = sprintf('$t$ = %1.1f h', t_/3600);
                set(l1, 'YData', z_, 'XData', thetaLine);
                data{k_, 1} = t_;
                data{k_, 2} = z_;
                data{k_, 3} = thetaLine;
                pause(rr);
            end
        end
        function [f, bulkTemp] = plotBulkVolumetricTemp(obj, dimensions)
            % plots the bulk volumetric temperature in the particle domain
            if nargin < 2, dimensions = 0; end
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            bulkTemp = zeros(length(obj.Fo), 1);
            t = zeros(size(obj.Fo));            
            z_ = obj.theta{2, 2};
            r_ = obj.theta{2, 3};
            theta_  = obj.theta{2, 1};
            fz = simpsonIntegrator(obj, z_);
            fr = simpsonIntegrator(obj, r_);
            Ir = 2*pi*(r_.*theta_)*fr';
            bulkTemp(1) = Ir'*fz'./(pi*obj.bp^2*max(z_));
            if dimensions, bulkTemp(1) = theta2T(obj, bulkTemp(1)); end
            for ii = 2:length(obj.Fo)
                ii_ = mod(ii-1, obj.ls) + 1;
                t(ii) = obj.Fo2t(obj.Fo(ii), 1);            
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                z_ = obj.theta{ii_, 2};
                r_ = obj.theta{ii_, 3};
                theta_  = obj.theta{ii_, 1};  
                fz = simpsonIntegrator(obj, z_);
                fr = simpsonIntegrator(obj, r_);
                Ir = 2*pi*(r_.*theta_)*fr';
                bulkTemp(ii) = Ir'*fz'./(pi*obj.bp^2*max(z_));
                if dimensions, bulkTemp(ii) = theta2T(obj, bulkTemp(ii)); end
            end
            f = figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.3], 'Visible', 'on');
            plot(t/3600, bulkTemp, '.-k');
            title('Bulk Volumetric Temperature', ...
                   'interpreter', 'latex', 'FontSize', 14);
            if dimensions
                ylabel('$T$ ($^\circ$C)', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.theta2T(obj.thetaA), max(obj.T0(:))]);
            else
                ylabel('$\theta$ ($^\circ$C)', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([0, 1]);
            end
            xlabel('$t$ (h)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'TickLabelInterpreter', 'latex')
            set(gcf, 'Color', [1 1 1])
        end
        function [f, t, wallTemps] = plotWallTemp(obj, prototype)
            % plots the transient temperature of each wall layer
            % load data
            N = size(obj.meshedWallInsulation, 1); % number of layers
            r_ = NaN*ones(N, 1);
            for i = 1:N, r_(i) = obj.meshedWallInsulation{i, 2}(2); end
            thetaW_ = NaN*ones(length(obj.Fo), N);
            t = NaN*ones(length(obj.Fo), 1);
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t(ii) = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t(ii) = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii == 2
                    loadTheta(obj, 1+obj.ls-1);
                end
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                thetaW_(ii, :) = obj.theta{ii_, 7};
            end
            wallTemps = obj.theta2T(thetaW_);
            % plot transient temperature
            f = figure('Units', 'normalized', 'color', 'white', ...
                                          'Position', [0 0 0.5 0.4]);
            colormap(bone);
            labels = cell(N, 1);
            for i = 1:N
                plot(t/3600, obj.theta2T(thetaW_(:, i)), 'LineWidth', 1); hold on;
                labels{i} = strcat(num2str(r_(i)), ' m');
            end
            legend(labels, ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Wall Temperature ($^\circ$C)', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14)
            xlim([0, max(t)/3600]);
        end
        function [f, t, wallFlux] = plotWallFlux(obj, prototype)
            % plots the transient temperature of each wall layer
            % load data
            N = size(obj.meshedWallInsulation, 1); % number of layers
            r_ = NaN*ones(N, 1);
            for i = 1:N, r_(i) = obj.meshedWallInsulation{i, 2}(2); end
            thetaW_ = NaN*ones(length(obj.Fo), N);
            wallFlux = NaN*ones(length(obj.Fo), N+1);
            t = NaN*ones(length(obj.Fo), 1);
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t(ii) = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t(ii) = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii == 2
                    loadTheta(obj, 1+obj.ls-1);
                end
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                thetaW_(ii, :) = obj.theta{ii_, 7};
                wallFlux(ii, 1) = obj.theta{ii_, 11};
            end
            wallTemps = obj.theta2T(thetaW_);
            for n_ = 2:N+1
                % compute fluxes according to RC circuit                
                if n_ == N+1
                    wallFlux(:, n_) = (wallTemps(:, n_-1) - obj.Tinf)./ ...
                                              obj.Rwall{n_, 2};
                else
                    wallFlux(:, n_) = (wallTemps(:, n_-1) - wallTemps(:, n_))./ ...
                                              obj.Rwall{n_, 2};
                end
            end
            % plot transient temperature
            f = figure('Units', 'normalized', 'color', 'white', ...
                                          'Position', [0 0 0.5 0.4]);
            colormap(bone);
            labels = cell(N, 1);
            for i = 1:N
                plot(t/3600, wallFlux(:, i), 'LineWidth', 1); hold on;
                labels{i} = strcat(num2str(r_(i)), ' m');
            end
            legend(labels, ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Wall Flux ($W/m^2$)', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14)
            xlim([0, max(t)/3600]);
        end       
        function [f, t, baseTemps] = plotBaseTemp(obj, prototype)
            % plots the transient temperature of each wall layer
            % load data
            N = size(obj.meshedBaseInsulation, 1); % number of layers
            z_ = NaN*ones(N, 1);
            for i = 1:N, z_(i) = obj.meshedBaseInsulation{i, 2}(2); end
            thetaB_ = NaN*ones(length(obj.Fo), N);
            t = NaN*ones(length(obj.Fo), 1);
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t(ii) = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t(ii) = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii == 2
                    loadTheta(obj, 1+obj.ls-1);
                end
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                thetaB_(ii, :) = obj.theta{ii_, 8};
            end
            baseTemps = obj.theta2T(thetaB_);
            % plot transient temperature
            f = figure('Units', 'normalized', 'color', 'white', ...
                                          'Position', [0 0 0.5 0.4]);
            colormap(bone);
            labels = cell(N, 1);
            for i = 1:N
                plot(t/3600, obj.theta2T(thetaB_(:, i)), 'LineWidth', 1); hold on;
                labels{i} = strcat(num2str(z_(i)), ' m');
            end
            legend(labels, ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Base Temperature ($^\circ$C)', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14)
            xlim([0, max(t)/3600]);
        end
        function [f, t, baseFlux] = plotBaseFlux(obj, prototype)
            % plots the transient temperature of each wall layer
            % load data
            N = size(obj.meshedBaseInsulation, 1); % number of layers
            z_ = NaN*ones(N, 1);
            for i = 1:N, z_(i) = obj.meshedBaseInsulation{i, 2}(2); end
            thetaB_ = NaN*ones(length(obj.Fo), N);
            baseFlux = NaN*ones(length(obj.Fo), N+1);
            t = NaN*ones(length(obj.Fo), 1);
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t(ii) = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t(ii) = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii == 2
                    loadTheta(obj, 1+obj.ls-1);
                end
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                thetaB_(ii, :) = obj.theta{ii_, 8};
                baseFlux(ii, 1) = obj.theta{ii_, 12};
            end
            baseTemps = obj.theta2T(thetaB_);
            for n_ = 2:N+1
                % compute fluxes according to RC circuit                
                if n_ == N+1
                    baseFlux(:, n_) = (baseTemps(:, n_-1) - obj.Tinf)./ ...
                                              obj.Rbase{n_, 2};
                else
                    baseFlux(:, n_) = (baseTemps(:, n_-1) - baseTemps(:, n_))./ ...
                                              obj.Rbase{n_, 2};
                end
            end
            % plot transient temperature
            f = figure('Units', 'normalized', 'color', 'white', ...
                                          'Position', [0 0 0.5 0.4]);
            colormap(bone);
            labels = cell(N, 1);
            for i = 1:N
                plot(t/3600, baseFlux(:, i), 'LineWidth', 1); hold on;
                labels{i} = strcat(num2str(z_(i)), ' m');
            end
            legend(labels, ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Base Flux ($W/m^2$)', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14)
            xlim([0, max(t)/3600]);
        end
        function [f, t, waTemps] = plotTopWallTemp(obj, prototype)
            % plots the transient temperature of each wall layer
            % load data
            N = size(obj.wallInsulation, 1); % number of layers
            r_ = NaN*ones(N, 1);
            for i = 1:N, r_(i) = obj.wallInsulation{i, 2}(2); end
            thetaW_ = NaN*ones(length(obj.Fo), N);
            t = NaN*ones(length(obj.Fo), 1);
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t(ii) = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t(ii) = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii == 2
                    loadTheta(obj, 1+obj.ls-1);
                end
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                thetaW_(ii, :) = obj.theta{ii_, 10};
            end
            waTemps = obj.theta2T(thetaW_);
            % plot transient temperature
            f = figure('Units', 'normalized', 'color', 'white', ...
                                          'Position', [0 0 0.5 0.4]);
            colormap(bone);
            labels = cell(N, 1);
            for i = 1:N
                plot(t/3600, obj.theta2T(thetaW_(:, i)), 'LineWidth', 1); hold on;
                labels{i} = strcat(num2str(r_(i)), ' m');
            end
            legend(labels, ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Top Wall Temperature ($^\circ$C)', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14)
            xlim([0, max(t)/3600]);
        end
        function [f, t, waFlux] = plotWallRadiation(obj, prototype)
            % plots the transient temperature of each wall layer
            % load data
            N = size(obj.wallInsulation, 1); % number of layers
            r_ = NaN*ones(N, 1);
            for i = 1:N, r_(i) = obj.wallInsulation{i, 2}(2); end
            thetaW_ = NaN*ones(length(obj.Fo), N);
            waFlux = NaN*ones(length(obj.Fo), N+1);
            t = NaN*ones(length(obj.Fo), 1);
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t(ii) = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t(ii) = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii == 2
                    loadTheta(obj, 1+obj.ls-1);
                end
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                thetaW_(ii, :) = obj.theta{ii_, 10};
                waFlux(ii, 1) = obj.theta{ii_, 18};
            end
            % plot transient temperature
            f = figure('Units', 'normalized', 'color', 'white', ...
                                          'Position', [0 0 0.5 0.4]);
            colormap(bone);
            labels = cell(N, 1);
            for i = 1:N
                plot(t/3600, waFlux(:, i), 'LineWidth', 1); hold on;
                labels{i} = strcat(num2str(r_(i)), ' m');
            end
            legend(labels, ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Wall Radiation ($W$)', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14)
            xlim([0, max(t)/3600]);
        end
        function [f, t, roofTemps] = plotRoofTemp(obj, prototype)
            % plots the transient temperature of each wall layer
            % load data
            N = size(obj.roofInsulation, 1); % number of layers
            t_ = NaN*ones(N, 1);
            for i = 1:N, t_(i) = obj.roofInsulation{i, 2}(2); end
            thetaR_ = NaN*ones(length(obj.Fo), N);
            t = NaN*ones(length(obj.Fo), 1);
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t(ii) = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t(ii) = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii == 2
                    loadTheta(obj, 1+obj.ls-1);
                end
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                thetaR_(ii, :) = obj.theta{ii_, 9};
            end
            roofTemps = obj.theta2T(thetaR_);
            % plot transient temperature
            f = figure('Units', 'normalized', 'color', 'white', ...
                                          'Position', [0 0 0.5 0.4]);
            colormap(bone);
            labels = cell(N, 1);
            for i = 1:N
                plot(t/3600, obj.theta2T(thetaR_(:, i)), 'LineWidth', 1); hold on;
                labels{i} = strcat(num2str(t_(i)), ' m');
            end
            legend(labels, ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Roof Layer Temperatures ($^\circ$C)', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14)
            xlim([0, max(t)/3600]);
        end
        function [f, t, roofFlux] = plotRoofRadiation(obj, prototype)
            % plots the transient temperature of each wall layer
            % load data
            N = size(obj.roofInsulation, 1); % number of layers
            t_ = NaN*ones(N, 1);
            for i = 1:N, t_(i) = obj.roofInsulation{i, 2}(2); end
            thetaR_ = NaN*ones(length(obj.Fo), N);
            roofFlux = NaN*ones(length(obj.Fo), N+1);
            t = NaN*ones(length(obj.Fo), 1);
            for ii = 2:length(obj.Fo)                                      
                ii_ = mod(ii-1, obj.ls) + 1;
                if prototype
                    t(ii) = obj.Fo2t(obj.Fo(ii), 1);
                else
                    t(ii) = obj.Fo2t(obj.Fo(ii));
                end
                % update domain according to mass accounting
                if ii == 2
                    loadTheta(obj, 1+obj.ls-1);
                end
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end
                thetaR_(ii, :) = obj.theta{ii_, 9};
                roofFlux(ii, 1) = obj.theta{ii_, 17};
            end
            % plot transient temperature
            f = figure('Units', 'normalized', 'color', 'white', ...
                                          'Position', [0 0 0.5 0.4]);
            colormap(bone);
            labels = cell(N, 1);
            for i = 1:N
                plot(t/3600, roofFlux(:, i), 'LineWidth', 1); hold on;
                labels{i} = strcat(num2str(t_(i)), ' m');
            end
            legend(labels, ...
                'interpreter', 'latex', 'FontSize', 16, ...
                'Location', 'northeast', 'NumColumns', 4)
            legend('boxoff')
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('Roof Radiation ($W$)', ... 
                   'interpreter', 'latex', 'FontSize', 16);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], 'OuterPosition', [0.01 0.09 0.9 0.9], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14)
            xlim([0, max(t)/3600]);
        end
        function f = plotOutletTempBulk(obj, dimensions)
            % plot outlet temperature
            f = figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.3], 'Visible', 'on');
            computeThetaO(obj);
            if dimensions
                plot(obj.thetaOB(:, 1)/3600, ...
                    obj.theta2T(obj.thetaOB(:, 2)), '.-k');
                title('Outlet Bulk Temperature: $\frac{2}{a^2}\int_0^aT_c\overline{r}d\overline{r}$', ...
                       'interpreter', 'latex', 'FontSize', 14);
                ylabel('$T_o$ ($^\circ$C)', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.theta2T(obj.thetaA), max(obj.T0(:))]);
                hold on
            else
                plot(obj.thetaOB(:, 1)/3600, obj.thetaOB(:, 2), '.-k');
                title('Outlet Bulk Temperature: $\frac{2}{a^2}\int_0^a\theta_c\overline{r}d\overline{r}$', ...
                   'interpreter', 'latex', 'FontSize', 14);           
                ylabel('$\theta_o$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.thetaA, max(obj.T2theta(obj.T0(:)))]);
                hold on
            end
            xlabel('$t$ (h)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'TickLabelInterpreter', 'latex')
            set(gcf, 'Color', [1 1 1])
            xlim([0, obj.thetaOB(end, 1)/3600]);
        end
        function f = compareOutletTempBulk(obj, f1, f2, f3, dimensions)
            % plot outlet temperature
            f = figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.3], 'Visible', 'on');
            obj.thetaFolder = f1;
            theta1_ = computeThetaO(obj);
            obj.thetaFolder = f2;
            theta2_ = computeThetaO(obj);
            obj.thetaFolder = f3;
            theta3_ = computeThetaO(obj);
            if dimensions
                plot(theta1_(:, 1)/3600, ...
                    obj.theta2T(theta1_(:, 2)), '-k', 'LineWidth', 1.5);
                title('Outlet Bulk Temperature: $\frac{2}{a^2}\int_0^aT_c\overline{r}d\overline{r}$', ...
                       'interpreter', 'latex', 'FontSize', 14);
                ylabel('$T$ ($^\circ$C)', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.theta2T(obj.thetaA), 810]);
                hold on
                plot(theta2_(:, 1)/3600, ...
                    obj.theta2T(theta2_(:, 2)), '-b', 'LineWidth', 1.5);
                plot(theta3_(:, 1)/3600, ...
                    obj.theta2T(theta3_(:, 2)), '-r', 'LineWidth', 1.5);
                fill([theta2_(:, 1)/3600; flipud(theta2_(:, 1)/3600)], ...
                    [obj.theta2T(theta2_(:, 2)); ...
                    flipud(obj.theta2T(theta3_(:, 2)))], ...
                    'g', 'FaceAlpha', 0.5);             
            else
                plot(obj.thetaOB(:, 1)/3600, obj.thetaOB(:, 2), '.-k');
                title('Outlet Bulk Temperature: $\frac{2}{a^2}\int_0^a\theta_c\overline{r}d\overline{r}$', ...
                   'interpreter', 'latex', 'FontSize', 14);           
                ylabel('$\theta_o$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.thetaA, max(obj.T2theta(obj.T0(:)))]);
                hold on
            end
            xlabel('$t$ (h)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'TickLabelInterpreter', 'latex')
            set(gcf, 'Color', [1 1 1])
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 14);
            xlim([16, 24]);
        end
        function [f, qw, qb, qt, qTot, energyIn, energyLoss, percentLoss] = plotBinLosses(obj, qwInc, qbInc, qtInc, qTotInc)
            % plots heat loss at discrete points in the composite wall
            % initialize domain ,temperature, and velocity distributions
            computeThetaO(obj);
            mtot = obj.QChp*obj.rhopPack*6*3600;
            qAve = -mtot*obj.cpp*(obj.theta2T(obj.thetaOB(end)) - 800)/(24*3600*1000);
            qTot = zeros(length(obj.Fo)-1, 1);
            qw = zeros(length(obj.Fo)-1, 1);
            qb = zeros(length(obj.Fo)-1, 1);
            qt = zeros(length(obj.Fo)-1, 1);
            t = zeros(length(obj.Fo)-1, 1);
            energyLoss = zeros(length(obj.Fo)-1, 1);
            energyIn = zeros(length(obj.Fo)-1, 1);
            percentLoss = zeros(length(obj.Fo)-1, 1);
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            for i = 2:length(obj.Fo)                                      
                i_ = mod(i-1, obj.ls) + 1;
                if i_ == 1 && i ~= length(obj.Fo)
                    loadTheta(obj, i+obj.ls-1);
                end
                t(i-1) = obj.Fo2t(obj.Fo(i), 1);
                if qwInc, qw(i-1) = obj.theta{i_, 11}; 
                else, qw(i-1) = NaN; end
                if qbInc, qb(i-1) = obj.theta{i_, 12};
                else, qb(i-1) = NaN; end
                if qtInc, qt(i-1) = obj.theta{i_, 13};
                else, qt = NaN; end
                if qTotInc, qTot(i-1) = obj.theta{i_, 14}; 
                else, qTot = NaN; end
                if i == 2
                    energyLoss(i-1) = 0;
                elseif i == 3
                    energyLoss(i-1) = trapz(t(1:i-1), qTot(1:i-1));
                else
                    ft = simpsonIntegrator(obj, t(1:i-1));
                    energyLoss(i-1) = ft*qTot(1:i-1);
                end
                if i == 2
                    energyIn(i-1) = obj.theta{i_, 15};
                else
                    energyIn(i-1) = energyIn(i-2) + obj.theta{i_, 15};
                end
                percentLoss(i-1) = 100*energyLoss(i-1)/energyIn(i-1);
            end
%             ft = simpsonIntegrator(obj, t);
%             energyLossTot = ft*qTot;           
%             energyInTot = mtot*obj.cpp*(obj.T0 - obj.Tinf)/1000;
%             percentLossTot = 100*energyLossTot/energyInTot;
            % plot total heat loss at each time step
            f = figure('Units', 'normalized', 'color', 'white', ...
                                     'Position', [0 0 0.5 0.3]); hold on;
            hold on;
            labels = cell(1, 4);
            if qTotInc 
                plot(t/3600, qTot, '-k'); 
                labels{1, 1} = 'Total Heat Loss'; 
            else
                labels{1, 1} = []; 
            end
            if qwInc 
                plot(t/3600, qw, '--b');
                labels{1, 2} = 'Wall Heat Loss'; 
            else
                labels{1, 2} = []; 
            end
            if qbInc
                plot(t/3600, qb, '--r'); 
                labels{1, 3} = 'Base Heat Loss'; 
            else
                labels{1, 3} = []; 
            end
            if qtInc
                plot(t/3600, qt, '--m'); 
                labels{1, 4} = 'Top Surface Heat Loss'; 
            else
                labels{1, 4} = []; 
            end 
%             plot(t/3600, qAve*ones(length(t), 1), '-r');
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$Q$ (kW)', 'interpreter', 'latex', ...
                                                      'FontSize', 14);
            legend(labels, 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 12);            
        end
        function f = plotOutletTempCL(obj, dimensions)
            % plot outlet temperature
            f = figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.3], 'Visible', 'on');  
            computeThetaO(obj);
            if dimensions
                plot(obj.thetaO(:, 1)/3600, ...
                    obj.theta2T(obj.thetaO(:, 2)), '.-k');
                title('Outlet Centerline Temperature', ...
                       'interpreter', 'latex', 'FontSize', 14);
                ylabel('$T_o$ ($^\circ$C)', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.theta2T(obj.thetaA), max(obj.T0(:))]);
                hold on
            else
                plot(obj.thetaO(:, 1)/3600, obj.thetaO(:, 2), '.-k');
                title('Outlet Centerline Temperature', ...
                   'interpreter', 'latex', 'FontSize', 14);           
                ylabel('$\theta_o$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.thetaA, max(obj.T2theta(obj.T0(:)))]);
                hold on
            end
            xlabel('$t$ (h)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'TickLabelInterpreter', 'latex')
            set(gcf, 'Color', [1 1 1])
            xlim([0, obj.thetaO(end, 1)/3600]);
        end
        function plotQWall(obj)
            % plots heat loss at discrete points in the composite wall
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            M = length(obj.wallInsulation(:, 1)) + 1;
            fz = simpsonIntegrator(obj, obj.zbarW);
            r_ = unique([obj.wallInsulation{:, 2}]);
            qw = cell(M, 1); t = [];
            for i = 2:length(obj.Fo)                                      
                i_ = mod(i-1, obj.ls) + 1;
                if i_ == 1 && i ~= length(obj.Fo)
                    loadTheta(obj, i+obj.ls-1);
                end
                t = [t, obj.Fo2t(obj.Fo(i), 1)];
                q_ = obj.theta{i_, 9};
                for ii = 1:M                   
                    qw{ii} = [qw{ii}, 2*pi*r_(ii)*(fz*q_{ii})*obj.Hp^2]; 
                end             
            end
            lstr = cell(1, M);
            for ii = 1:M, lstr{ii} = sprintf('r = %1.2f', r_(ii)); end
            % plot total heat transfer for each segment
            figure('Units', 'normalized', 'color', 'white', ...
                                     'Position', [0 0 0.5 0.3]); hold on;
            for ii = 1:M 
                lstr{ii} = sprintf('r = %1.2f', r_(ii));
                plot(t, qw{ii});
            end
            xlabel('$t (s)$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$Q_{wall}$ (W)', 'interpreter', 'latex', ...
                                                      'FontSize', 14);
            legend(lstr);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 12) 
            
        end
        function plotAirTemp(obj)
            % plots heat loss at discrete points in the composite wall
            % initialize domain ,temperature, and velocity distributions
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            t = zeros(length(obj.Fo), 1);
            Ta_ = zeros(length(obj.Fo), 1);
            for i = 2:length(obj.Fo)                                      
                i_ = mod(i-1, obj.ls) + 1;
                if i_ == 1 && i ~= length(obj.Fo)
                    loadTheta(obj, i+obj.ls-1);
                end
                t(i) = obj.Fo2t(obj.Fo(i), 1);
                Ta_(i) = obj.theta2T(obj.theta{i_, 12});             
            end
            % plot total heat transfer for each segment
            figure('Units', 'normalized', 'color', 'white', ...
                                     'Position', [0 0 0.5 0.3]); hold on;
            plot(t(2:end-1)/3600, Ta_(2:end-1), '-k');
            xlabel('$t (h)$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$T_a$ ($^\circ C$)', 'interpreter', 'latex', ...
                                                      'FontSize', 14);
            set(gca, 'box', 'off', 'TickDir', 'both', ...
                'TickLength', [0.01, 0.025], ...
                'TickLabelInterpreter', 'latex', 'FontSize', 12) 
            
        end
        function plotTempPoint(obj, z_, r_, dimensions, exp_cmp)
            % plots temperature for specified zbar and rbar
            figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.3], 'Visible', 'on');
            title(sprintf('Temperature at $z/H$ = %1.2f and $r/H$ = %1.2f', ...
                       z_, r_), 'interpreter', 'latex', 'FontSize', 14);
            hold on
            theta_ = NaN*ones(length(obj.Fo));
            [~, i_] = min(abs(z_ - obj.z));
            [~, j_] = min(abs(r_ - obj.r));
            for ii = 1:length(obj.Fo)
                ii_ = mod(ii-1, obj.ls) + 1;
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);               
                end
                thetak_ = thetak(obj, ii_);
                theta_(ii) = thetak_(i_, j_);
            end             
            if dimensions
                plot(obj.Fo2t(obj.Fo), obj.theta2T(theta_), '.-k');
                ylabel('$T$ ($^\circ$C)', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.theta2T(obj.thetaA), obj.T0]);
            else         
                plot(obj.Fo2t(obj.Fo), theta_, '.-k');
                ylabel('$\theta$', 'interpreter', 'latex', ...
                    'FontSize', 14);
                ylim([obj.thetaA, obj.T2theta(obj.T0)]);
            end
            xlabel('$t$ (s)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'TickLabelInterpreter', 'latex')
            set(gcf, 'Color', [1 1 1])
            if exp_cmp               
                if isempty(obj.thetaExp)
                    warning('no experimental data found for comparison') 
                else
                    thetaExp_ = NaN*ones(length(obj.Fo));
                    for ii = 1:length(obj.Fo)
                        thetakExp_ = thetakExp(obj, ii);
                        thetaExp_(ii) = thetakExp_(i_, j_);
                    end 
                    if dimensions
                        plot(obj.Fo2t(obj.Fo), obj.theta2T(thetaExp_), 'xr');
                    else
                        plot(obj.Fo2t(obj.Fo), thetaExp_, ...
                            'xr');
                    end
                end                
            end
        end
        function plotFitResults(obj, result)
            % plots experimental comparison for fitting results contained
            % in the result input provided
            for i = 1:length(result.thetaFit)
                figure('Units', 'normalized', ...
                        'Position', [0 0 0.4 0.3], 'Visible', 'on');
                title(sprintf('Temperature at $z/H$ = %1.2f and $r/H$ = %1.2f', ...
                   result.thetaFit{i}.z, result.thetaFit{i}.r), ...
                   'interpreter', 'latex', 'FontSize', 14);
                hold on
                plot(obj.Fo2t(obj.Fo),  ...
                    obj.theta2T(result.thetaFit{i}.theta), '.k');
                [~, ie] = min(abs(obj.zStore(i) - obj.zbarExp) ...
                    + abs(obj.rStore(i) - obj.rbarExp));
                plot(obj.Fo2t(obj.FoExp),  ...
                    obj.theta2T(result.thetaExp{ie}.theta), 'xr');
                legend('Simulation', 'Experiment', ...
                    'interpreter', 'latex', 'FontSize', 14)
            end
        end
        function [fzri, pzri1, pzri2] = plotZRTemp(obj, temp, z, r, tempW, ...
                                                     zw_, rw_, ShowPlot)
            % plots ZR plane for one time instance
            [R, Z] = meshgrid(r, z); 
            [RW, ZW] = meshgrid(rw_, zw_);
            % plot contour
            fzri = figure('Units', 'normalized', ...
                'Position', [0 0 obj.b 1], 'Visible', 'off');
            pzri1 = surf(R, Z, temp, 'EdgeColor', 'interp', ...
                 'FaceColor', 'interp'); hold on
            pzri2 = surf(RW, ZW, tempW);
            xlim([0, max(RW(:))])
            ylim([min(Z(:)), max(ZW(:))])
            pbaspect([obj.b, 1, 1]);  % figure sized proportional to aspect ratio
            view(0, 90);
            caxis([obj.T2theta(obj.Tinf), 1]);
            colormap(jet);
            cb = colorbar;
            cb.Ruler.MinorTick = 'on';
            set(pzri1, 'linestyle', 'none');
            set(pzri2, 'linestyle', 'none');
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
            if nargin > 2 && ShowPlot
                fzri.Visible = 'on';
            end
        end
        function [fzri, pzri1] = plotZRTempNoWall(obj, temp, z, r, ShowPlot)
            % plots ZR plane for one time instance
            [R, Z] = meshgrid(r, z); 
            % plot contour
            fzri = figure('Units', 'normalized', ...
                'Position', [0 0 2*obj.b 1], 'Visible', 'off');
            pzri1 = surf(R, Z, temp, 'EdgeColor', 'interp', ...
                 'FaceColor', 'interp');
            xlim([min(R(:)), max(R(:))])
            ylim([min(Z(:)), max(Z(:))])
            pbaspect([2*obj.b, 1, 1]);  % figure sized proportional to aspect ratio
            view(0, 90);
            caxis([obj.T2theta(obj.Tinf), 1]);
            colormap(jet);
            cb = colorbar;
            cb.Ruler.MinorTick = 'on';
            set(pzri1, 'linestyle', 'none');
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
            if nargin > 2 && ShowPlot
                fzri.Visible = 'on';
            end
        end
        function [fzri, pzri] = plotZRCenterSensitivity(obj, s, z, r, ...
                                  ShowPlot)
            % plots ZR plane for one time instance
            [R_, Z_] = meshgrid(r, z); 
            [R, Z] = meshgrid(r, linspace(0, z(end), size(s, 1)));
            s_ = interp2(R, Z, s, R_, Z_);
            % plot contour
            fzri = figure('Units', 'normalized', ...
                'Position', [0 0 4*obj.a0 1], 'Visible', 'off');
            pzri = surf(R_, Z_, s_/max(s_(:)), 'EdgeColor', 'interp', ...
                'FaceColor', 'interp');         
            xlim([0, max(R_(:))])
            ylim([min(Z_(:)), max(Z_(:))])
            pbaspect([4*obj.a0, 1, 1]);  % figure sized proportional to aspect ratio
            view(0, 90);
            caxis([0, 1]);
            colormap(magma);
            cb = colorbar;
            cb.Ruler.MinorTick = 'on';
            set(pzri, 'linestyle', 'none');
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
            if nargin > 2 && ShowPlot
                fzri.Visible = 'on';
            end
        end 
        function [fzri, pzri] = plotZRConservation(obj, c, z, r, ShowPlot)
            % plots ZR plane for one time instance
            [R, Z] = meshgrid(r, z); 
            % plot contour
            fzri = figure('Units', 'normalized', ...
                'Position', [0 0 obj.b 1], 'Visible', 'off');
            pzri = surf(R, Z, c, 'EdgeColor', 'interp', ...
                                               'FaceColor', 'interp');         
            xlim([0, max(R(:))])
            ylim([min(Z(:)), max(Z(:))])
            pbaspect([obj.b, 1, 1]);  % figure sized proportional to aspect ratio
            view(0, 90);
            caxis([1e-6, 1]);
            colormap(flipud(magma));
            cb = colorbar;
            cb.Ruler.MinorTick = 'on';
            set(pzri, 'linestyle', 'none');
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
            if nargin > 2 && ShowPlot
                fzri.Visible = 'on';
            end
        end   
        function [fzri, pzri] = plotZRPsi(obj, temp, ShowPlot)
            % plots ZR plane for one time instance
            z_ = obj.zbar;
            z0_ = obj.zbar0;                        
            [R, Z] = meshgrid(obj.rbar, z_); 
            R = R(1:length(z0_), :);
            Z = Z(1:length(z0_), :);
            % plot contour
            fzri = figure('Units', 'normalized', ...
                'Position', [0 0 obj.b 1], 'Visible', 'off');
            pzri = surf(R, Z, temp);         
            xlim([0, max(R(:))])
            ylim([min(Z(:)), max(Z(:))])
            pbaspect([obj.b, 1, 1]);  % figure sized proportional to aspect ratio
            view(0, 90);
            caxis([obj.thetaA, 1]);
            colormap(flipud(autumn));
            cb = colorbar;
            cb.Ruler.MinorTick = 'on';
            set(pzri, 'linestyle', 'none');
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
            if nargin > 2 && ShowPlot
                fzri.Visible = 'on';
            end
        end  
        function captureDimensionalTemp(obj, k_)
            % plots the dimensional temperature profile over the entire
            % domain at the indicated time step, k
            temp = obj.theta2T(obj.thetak(k_));
            z_ = obj.theta{k_, 2};
            r_ = obj.theta{k_, 3};
            t = obj.Fo2t(obj.Fo(k_), 1);
            [R, Z] = meshgrid(r_, z_); 
            % plot contour
            fzri = figure('Units', 'normalized', ...
                'Position', [0 0 obj.b 1], 'Visible', 'off');
            pzri = surf(R, Z, temp); %, 'EdgeColor', 'interp', ...
%                 'FaceColor', 'interp');         
            xlim([0, max(R(:))])
            ylim([min(Z(:)), max(Z(:))])
            pbaspect([obj.b, 1, 1]);  % figure sized proportional to aspect ratio
            view(0, 90);
            caxis([700, 800]);
            colormap(jet);
            cb = colorbar;
            cb.Ruler.MinorTick = 'on';
            set(pzri, 'linestyle', 'none');
            ylabel(cb, '$T$ ($^\circ$C)', 'interpreter', 'latex', 'FontSize', 14);
            xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
            ylabel('$z/H$', 'interpreter', 'latex', 'FontSize', 14);
            title(sprintf('$t$ = %1.0f h', t/3600), 'interpreter', 'latex', ...
                'FontSize', 14);
            set(gca, 'TickLabelInterpreter', 'latex')
            set(gcf, 'Color', [1 1 1])
            if nargin > 2 && ShowPlot
                fzri.Visible = 'on';
            end      
            % save figure
            saveas(pzri, sprintf('tempDegC_%1.0d.png', k_));
        end
        function [fzri, pzri] = plotZRVel(obj, vel, ShowPlot)
            % plots ZR plane for one time instance 
            [R, Z] = meshgrid([obj.rhat, obj.rbar], ...
                [obj.zbar, (1 + obj.zhat)]);                   
            % plot contour
            fzri = figure('Units', 'normalized', ...
                'Position', [0 0 obj.b 1], 'Visible', 'off');
            pzri = surf(R, Z, vel);         
            xlim([0, max(R(:))])
            ylim([min(Z(:)), max(Z(:))])
            pbaspect([obj.b, 1, 1]);  % figure sized proportional to aspect ratio
            view(0, 90);
            caxis([0, 1]);
            colormap(flipud(winter));
            cb = colorbar;
            cb.Ruler.MinorTick = 'on';
            set(pzri, 'linestyle', 'none');
            ylabel(cb, '$|U|/U_\infty$', 'interpreter', 'latex', 'FontSize', 24);
            xlabel('$\overline{r}$', 'interpreter', 'latex', 'FontSize', 24);
            ylabel('$\overline{z}$', 'interpreter', 'latex', 'FontSize', 24);
            title(sprintf('Velocity Magnitude'), 'interpreter', ...
                'latex', 'FontSize', 14);
            if nargin > 2 && ShowPlot
                fzri.Visible = 'on';
            end
        end
        function plotMassAccounting(obj)
            % plot outlet temperature
            if isempty(obj.thryMass)
                computeTheoreticalMass(obj)
            end
            figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.3], 'Visible', 'on');            
            plot(obj.Fo2t(obj.thryMass(:, 1)), ...
                   obj.thryMass(:, 2), '.-k');
            title('Total Mass in Storage Bin', 'interpreter', 'latex', ...
                'FontSize', 14);
            ylabel('$m$ (kg)', 'interpreter', 'latex', 'FontSize', 14);
            ylim([0, max(obj.thryMass(:, 2))]);
            xlabel('$t$ (s)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'TickLabelInterpreter', 'latex')
            set(gcf, 'Color', [1 1 1])
            hold on
            if ~isempty(obj.simMass)
                plot(obj.Fo2t(obj.simMass(:, 1)), obj.simMass(:, 2), 'x-r');
                legend('Theoretical', 'Simulation', ...
                                  'interpreter', 'latex', 'FontSize', 14)
            end
        end
        function plotTopLoss(obj)
            % plots the radial temperature change for each vertical
            % cell over the simulation
            figure('Units', 'normalized', ...
                'Position', [0 0 0.4 0.3], 'Visible', 'on');    
            hold on
            plot(obj.Fo2t(obj.topLoss{1}), obj.topLoss{3}, '-k');
            title('Top Convective Losses', 'interpreter', 'latex', ...
                'FontSize', 14);
            ylabel('$\Delta T$ ($^\circ C$)', 'interpreter', 'latex', ...
                'FontSize', 14);
            xlabel('$t$ (s)', 'interpreter', 'latex', 'FontSize', 14);
            set(gca, 'TickLabelInterpreter', 'latex')
            set(gcf, 'Color', [1 1 1])
            hold on
        end            
        function snapContour(obj, snapK, HoldTime)
            if nargin < 3, HoldTime = 0; end
            % captures a sequence of still contour plots of temperature
            % distribution
            figure('Units', 'normalized', 'Position', [0 0 1*obj.b 0.5]);
            snap = NaN*ones(1, snapK);
            Fo_ = linspace(obj.t2Fo(100), obj.t2Fo(300), snapK);
            for i = 1:snapK
                [~, snap(i)] = min(abs(obj.Fo - Fo_(i)));
            end
            if HoldTime ~= 0
                [~, kstart] = min(abs(obj.t2Fo(HoldTime) - obj.Fo));
                loadTheta(obj, kstart-mod(kstart-1, obj.ls)+obj.ls);
            else
                kstart = 1; loadTheta(obj, obj.ls);
            end
            k_ = 1;
            for ii = kstart:length(obj.Fo)-1
                ii_ = mod(ii-1, obj.ls) + 1;
                t = obj.Fo2t(obj.Fo(ii)) - HoldTime;
                % update domain according to mass accounting
                if ii_ == 1 && ii ~= length(obj.Fo)
                    loadTheta(obj, ii+obj.ls-1);
                end                
                z_ = [obj.theta{ii_, 2}, ...
                                  max(obj.theta{ii_, 2}), 1 + 1.000001*obj.h];
                r_ = obj.theta{ii_, 3};
                theta_  = [obj.theta{ii_, 1}; obj.thetaA*zeros(2, length(r_))];
                if snap(k_) == ii
                    subplot(1, snapK, k_);
                    % plots ZR plane for one time instance
                    [R, Z] = meshgrid(r_, z_); 
                    % plot contour
                    pzri = surf(R, Z, obj.theta2T(theta_), ...
                        'EdgeColor', 'interp', 'FaceColor', 'interp');         
                    xlim([0, max(R(:))])
%                     ylim([min(Z(:)), max(Z(:))])
                    ylim([min(Z(:)), 0.6])
                    pbaspect([obj.b, 1, 1]);
                    view(0, 90);
                    colormap(jet);
                    cb = colorbar;
                    set(pzri, 'linestyle', 'none');
                    xlabel('$r/H$', 'interpreter', 'latex', 'FontSize', 14);
                    xtickangle(45);                    
                    set(gca, 'box', 'off', 'TickDir', 'both', ...
                        'TickLength', [0.01, 0.025], ...
                        'TickLabelInterpreter', 'latex', 'FontSize', 12); 
                    caxis([600, 790]);
                    title(sprintf('$t$ = %1.0f s', t), ...
                                'interpreter', 'latex', 'FontSize', 12);                     
                    if k_ == 1
                        ylabel('$z/H$', 'interpreter', 'latex', ...
                            'FontSize', 14);
%                         title(sprintf('$t$ = %1.0f s', 0), ...
%                                 'interpreter', 'latex', 'FontSize', 12);
                    else
                        set(gca, 'YTickLabel', []);
                        pos = get(gca, 'Position');
                        pos(1) = (0.9 - 0.015*k_)*pos(1);
                        set(gca, 'Position', pos);
                    end
                    if k_ == snapK
                        cb.Ruler.MinorTick = 'on';
                        ylabel(cb, '$T$ ($^\circ$C)', ...
                               'interpreter', 'latex', 'FontSize', 14);
                        set(cb, 'box', 'off', 'TickDir', 'both', ...
                        'TickLabelInterpreter', 'latex', 'FontSize', 12); 
                    else
                        set(cb, 'YTickLabel', []);
                    end
                    k_ = k_ + 1;
                    if k_ > length(snap), break; end
                end
            end
        end
        function plotTopChannel(obj)
            % plots zoomed-in contour of the top flow channel at specified
            % time
            figure('Units', 'normalized', 'Position', [0 0 obj.b 0.05]);
%             ii_ = mod(length(obj.Fo)-1, obj.ls) + 1;
            ii_ = 3;
            % update domain according to mass accounting
%             loadTheta(obj, length(obj.Fo));   
            loadTheta(obj, obj.ls); 
            z_ = obj.theta{ii_, 2}(end-15:end);
            r_ = obj.theta{ii_, 3};
            theta_  = obj.theta{ii_, 1}(end-15:end, :);
            % plots ZR plane for one time instance
            [R, Z] = meshgrid(r_, z_); 
            % plot contour
            pzri = surf(R, Z, obj.theta2T(theta_), ...
                'EdgeColor', 'interp', 'FaceColor', 'interp');         
            xlim([0, max(R(:))])
            ylim([min(Z(:)), max(Z(:))])
%             pbaspect([obj.b, 1, 1]);
            view(0, 90);
            colormap(jet);
            set(pzri, 'linestyle', 'none');
            caxis([500, 800]);                  
            set(gca, 'YTickLabel', [], 'XTickLabel', []);
        end
        function plotCenterChannel(obj)
            % plots zoomed-in contour of the top flow channel at specified
            % time
            figure('Units', 'normalized', 'Position', [0 0 0.5*obj.a0 obj.H]);
%             ii_ = mod(length(obj.Fo)-1, obj.ls) + 1;
            ii_ = 21;
            % update domain according to mass accounting
   
            loadTheta(obj, obj.ls); 
            z_ = obj.theta{ii_, 2};
            r_ = obj.theta{ii_, 3}(1:18);
            theta_  = obj.theta{ii_, 1}(:, 1:18);
            % plots ZR plane for one time instance
            [R, Z] = meshgrid(r_, z_); 
            % plot contour
            pzri = surf(R, Z, obj.theta2T(theta_), ...
                'EdgeColor', 'interp', 'FaceColor', 'interp');         
            xlim([0, max(R(:))])
            ylim([min(Z(:)), max(Z(:))])
%             pbaspect([obj.b, 1, 1]);
            view(0, 90);
            colormap(jet);
            set(pzri, 'linestyle', 'none');
            caxis([600, 790]);                  
            set(gca, 'YTickLabel', [], 'XTickLabel', []);
        end
        function saveData(obj)
            % saves simulation data to spreadsheet
            if length(obj.Fo) < obj.ls, loadTheta(obj, length(obj.Fo));  
            else, loadTheta(obj, obj.ls); end
            computeThetaO(obj);
            t = zeros(length(obj.Fo), 1);
            Tout = [0; obj.thetaOB(:, 2)];
            TvAve = zeros(length(obj.Fo), 1);
            qLossTot = zeros(length(obj.Fo), 1);            
            Ta_ = zeros(length(obj.Fo), 1);
            z_ = obj.theta{2, 2};
            r_ = obj.theta{2, 3};
            theta_  = obj.theta{2, 1};
            fz = simpsonIntegrator(obj, z_);
            fr = simpsonIntegrator(obj, r_);
            Ir = 2*pi*(r_.*theta_)*fr';
            TvAve(1) = obj.theta2T(Ir'*fz'./(pi*obj.bp^2*max(z_)));
            for i = 2:length(obj.Fo)                                      
                i_ = mod(i-1, obj.ls) + 1;
                if i_ == 1 && i ~= length(obj.Fo)
                    loadTheta(obj, i+obj.ls-1);
                end
                t(i) = obj.Fo2t(obj.Fo(i), 1);                
                qLossTot(i) = obj.theta{i_, 14};
                Ta_(i) = obj.theta2T(obj.theta{i_, 12});
                z_ = obj.theta{i_, 2};
                r_ = obj.theta{i_, 3};
                theta_  = obj.theta{i_, 1};  
                fz = simpsonIntegrator(obj, z_);
                fr = simpsonIntegrator(obj, r_);
                Ir = 2*pi*(r_.*theta_)*fr';
                TvAve(i) = obj.theta2T(Ir'*fz'./(pi*obj.bp^2*max(z_)));
            end
            % save data to spreadsheet
            tbl = table(t, obj.theta2T(Tout), TvAve, qLossTot, Ta_);
            writetable(tbl, 'FFSimTable.xlsx');                       
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % other
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        function theta_ = thetak(obj, k_)
            % returns theta matrix at time step k_
            theta_ = full(obj.theta{k_, 1});
        end       
        function saveTheta(obj, k_)
            % saves theta matrix at time step k_
            if k_ == length(obj.Fo)
                thetaString = sprintf('theta_%d.mat', floor(k_/obj.ls) + 1);
            else
                thetaString = sprintf('theta_%d.mat', floor(k_/obj.ls));
            end
            saveString = strcat(obj.thetaFolder, '\', thetaString);
            thetaSave = obj.theta;
            save(saveString, 'thetaSave');  
        end
        function saveThetaK(obj, k_)
            % saves theta matrix at time step k_
            if k_ == length(obj.Fo)
                thetaString = sprintf('thetaK_%d.mat', floor(k_/obj.ls) + 1);
            else
                thetaString = sprintf('thetaK_%d.mat', floor(k_/obj.ls));
            end
            saveString = strcat(obj.thetaFolder, '\', thetaString);
            thetaSave = obj.thetaK;
            save(saveString, 'thetaSave');  
        end
        function loadTheta(obj, k_)
            % loads theta matrix from stored location. The range of data
            % contained in the loaded matrix ends at time step k_
            if k_ == length(obj.Fo)
                thetaFile = sprintf('theta_%d.mat', floor(k_/obj.ls) + 1);
            else
                thetaFile = sprintf('theta_%d.mat', floor(k_/obj.ls));
            end
            thetaPath = strcat(obj.thetaFolder, '\', thetaFile);
            load(thetaPath, 'thetaSave');
            obj.theta = thetaSave;
        end 
        function loadThetaK(obj, k_)
            % loads theta matrix from stored location. The range of data
            % contained in the loaded matrix ends at time step k_
            if k_ == length(obj.Fo)
                thetaFile = sprintf('thetaK_%d.mat', floor(k_/obj.ls) + 1);
            else
                thetaFile = sprintf('thetaK_%d.mat', floor(k_/obj.ls));
            end
            thetaPath = strcat(obj.thetaFolder, '\', thetaFile);
            load(thetaPath, 'thetaSave');
            obj.thetaK = thetaSave;
        end 
        function tb = bulkTempR(obj, t, r_)
            % computes the bulk temperature for an r-dimensional array
            fr = simpsonIntegrator(obj, r_);
            tb = (t.*r_)*fr'/(r_*fr');
        end
        function tb = bulkTempZ(~, t)
            % computes the bulk temperature for an r-dimensional array
            tb = mean(t, 1);
        end
        function [qw, qb, qt, qTot] = computeBinHeatLoss(obj, thetaWall, ...
                                         drw, zw, thetaBase, dzb, rb, ...
                                         thetaTop, dzt, rt)
            % computes the heat loss rate at the base, wall, and top
            % surface at a particular time step with Fourier's Law
            qw = -obj.kp*(obj.T0 - obj.Tinf)*(thetaWall(:, 2) ...
                - thetaWall(:, 1))./drw;
            qb = obj.kp*(obj.T0 - obj.Tinf)*(thetaBase(2, :) ...
                - thetaBase(1, :))./dzb;
            qt = -obj.kp*(obj.T0 - obj.Tinf)*(thetaTop(2, :) ...
                - thetaTop(1, :))./dzt;
            % integrate to get total heat loss at each surface
            fzw = simpsonIntegrator(obj, zw);
            frb = simpsonIntegrator(obj, rb);
            frt = simpsonIntegrator(obj, rt);
            qw = (2*pi*obj.bp*obj.Hp^2*fzw*qw)/1000;       % kW
            qb = (2*pi*obj.Hp^2*frb*(qb.*rb)')/1000;       % kW
            qt = (2*pi*obj.Hp^2*frt*(qt.*rt)')/1000;       % kW
            qTot = abs(qw) + abs(qb) + abs(qt);                                                         
        end
        function c = variableC(~, t)
            % computes the temperature-dependent specific heat for every
            % value in t. Outputs specific heat matrix of same dim as t. t
            % must be in Kelvin.
            a_ = 148.2; b_ = 0.3093; c = a_*t.^b_;
        end
        function w = simpsonWeights(~, p)
            % generates a set of weights from an array of points for the
            % 3-point quadrature rule used in the Simpson Integrator
            w = NaN*ones(size(p)); n = size(p, 1);
            % build weighting matrix neglecting first row
            a_ = p(:, 1); b_ = p(:, 2); c_ = p(:, 3);
            det = 1./(a_.^2.*(c_-b_) + b_.^2.*(a_-c_) + c_.^2.*(b_-a_));
            w(:, 1) = det.*((b_ - c_).*(b_.^3 - c_.^3)/3 + ...
                            (b_ - c_).*(b_.^2.*c_ - b_.*c_.^2) + ...
                            (c_.^2 - b_.^2).*(b_.^2 - c_.^2)/2);
            w(:, 2) = det.*((a_.^2 - c_.^2).*(b_.^2 - c_.^2)/2 + ...
                            (b_ - c_).*(a_.*c_.^2 - a_.^2.*c_) + ...
                            (c_ - a_).*(b_.^3 - c_.^3)/3);
            w(:, 3) = det.*((a_ - b_).*(b_.^3 - c_.^3)/3 + ...
                            (b_ - c_).*(a_.^2.*b_ - a_.*b_.^2) + ...
                            (b_.^2 - a_.^2).*(b_.^2 - c_.^2)/2);
            w = spdiags(w, [0, 1, 2], n, n+2);
            % add first row with reversed quadrature
            w = [zeros(1, n+2); w]; 
            w(1, 1) = w(2, 3); w(1, 2) = w(2, 2); w(1, 3) = w(2, 1); 
        end
        function x = simpsonIntegrator(obj, v)
            % generates a vector that can be used to integrate a matrix or
            % vector along the dimension with a spatial mesh given by v
            n = length(v); v = reshape(v, n, 1);
            p_ = [v(1:end-2), v(2:end-1), v(3:end)];
            w = simpsonWeights(obj, p_);  
            x = sum(w);
        end
        function [x, dv] = simpsonIntegratorEq(~, v)
            % generates a vector that can be used to integrate a matrix or
            % vector along the dimension with a spatial mesh given by v
            dv = abs(v(2) - v(1));
            x = NaN*ones(size(v));
            x(1) = dv/3; x(2) = 15*dv/12; x(3) = 11*dv/12;
            x(4:end-2) = dv; x(end-1) = 13*dv/12; x(end) = 5*dv/12;
        end
        function [x, dx] = nodeGen(~, xlim, n)
            % generates a set of n chebyshev nodes spaced between xlim(1)
            % and xlim(2)
%             r_ = (xlim(2) - xlim(1))/2; theta_ = linspace(pi, 0, n);
%             x = xlim(1) + r_*(1 + cos(theta_));
%             dx = (eye(n) - diag(ones(1, n-1), -1))*x'; dx = dx(2:end);
              x = linspace(xlim(1), xlim(2), n);
              dx = x(2) - x(1);
        end
        function D = diffR(~, n, m, dr_)
            % creates a differential operator for first derivative in r 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;
            % set central differencing for center points
            D = spdiags([ones(nm, 1), -ones(nm, 1)], [1, -1], nm, nm);
            D = D./(2*dr_);
            % delete boundary entries
            D(m:m:end, :) = 0;      % r = 0
            D(m+1:m:end, :) = 0;    % r = b            
            % set forward differencing for left boundary
            D(1:m*(nm+1):end) = -1/dr_;
            D(nm+1:m*(nm+1):end) = 1/dr_;           
            % set backward differencing for right boundary
            D((m-1)*(nm+1)+1:m*(nm+1):end) = 1/dr_;
            D((m-2)*(nm+1)+2:m*(nm+1):end) = -1/dr_; 
        end   
        function D = diff2R(~, n, m, dr_)
            % creates a differential operator for first derivative in r 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;            
            % set central differencing for center points
            D = spdiags([-2*ones(nm, 1), ones(nm, 1), ones(nm, 1)], ...
                             [0, 1, -1], nm, nm);
            D = D./dr_^2;            
            % delete boundary entries
            D(m:m:end, :) = 0;      % r = 0
            D(m+1:m:end, :) = 0;    % r = b            
            % set forward differencing for left boundary
            D(1:m*(nm+1):end) = 1/dr_^2;
            D(nm+1:m*(nm+1):end) = -2/dr_^2;
            D(2*nm+1:m*(nm+1):end) = 1/dr_^2; 
            % set backward differencing for right boundary
            D((m-1)*(nm+1)+1:m*(nm+1):end) = 1/dr_^2;
            D((m-2)*(nm+1)+2:m*(nm+1):end) = -2/dr_^2;
            D((m-3)*(nm+1)+3:m*(nm+1):end) = 1/dr_^2; 
        end
        function D = diffZ(~, n, m, dz_)
            % creates a differential operator for first derivative in z 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;            
            % set central differencing for center points
            D = spdiags([ones(nm, 1), -ones(nm, 1)], [m, -m], nm, nm);
            D = D./(2*dz_);            
            % delete boundary entries
            D(1:m, :) = 0;          % z = 0
            D(end-m+1:end, :) = 0;    % z = 1            
            % set forward differencing for top boundary
            D(1:nm+1:m*(nm+1)) = -1/dz_;
            D(m*nm+1:nm+1:2*m*nm+1) = 1/dz_;            
            % set backward differencing for bottom boundary
            D(end-(m-1)*(nm+1):nm+1:end) = 1/dz_;
            D(end-2*m*(nm+1)+nm+1+m:nm+1:end-m*nm) = -1/dz_;               
        end
        function D = diff2Z(~, n, m, dz_)
            % creates a differential operator for first derivative in z 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;            
            % set central differencing for center points
            D = spdiags([-2*ones(nm, 1), ones(nm, 1), ones(nm, 1)], ...
                          [0, m, -m], nm, nm);
            D = D./dz_^2;            
            % delete boundary entries
            D(1:m, :) = 0;            % z = 0
            D(end-m+1:end, :) = 0;    % z = 1            
            % set forward differencing for top boundary
            D(1:nm+1:m*(nm+1)) = 1/dz_^2;
            D(m*nm+1:nm+1:2*m*nm+1) = -2/dz_^2;
            D(2*m*nm+1:nm+1:3*m*nm+1) = 1/dz_^2;   
            % set backward differencing for bottom boundary
            D(end-(m-1)*(nm+1):nm+1:end) = 1/dz_^2;
            D(end-2*m*(nm+1)+nm+1+m:nm+1:end-m*nm) = -2/dz_^2;
            D(end-3*m*(nm+1)+nm+1+2*m:nm+1:end-2*m*nm) = 1/dz_^2;
        end      
        % static analytic solutions         
        function t = Xt(~, Fo_, beta_, eta_)
            % transient component of analytic solutions
            t = exp(-(beta_^2 + eta_^2)*Fo_);
        end
    end
    methods(Static)
    end
end