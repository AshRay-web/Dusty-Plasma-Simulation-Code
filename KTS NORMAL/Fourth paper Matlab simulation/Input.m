%Prepared By NORMAL Roshan Chalise at 05/09/2012
% To study the sheath structure for an oblique magnetic field:
    Angle =45;
    Theta=(Angle*pi/180);                               % In terms of radian

% Basic parameters:
    B0=0.25;           %-300e-3                                          % Applied magnetic field(Tesla)
    vy=-1.0e+005;                                          % Value of Y-component velocity
    vz=-1.0e+005;                                           % Value of Z-component
    kB=1.38062e-23;                                     % Boltzmann constant(J/K)
    epsilon=8.85419e-12;                                % Permittivity of the medium(F/m)
    e=1.602192e-19;                                     % Electronic charge(C)
    Ze=e;                                               % Ion charge(C)
    Mi=1.672e-27;                                       % Mass of ion(Kg)
    Me=9.109e-31;                                       % Mass of electron(Kg)
    mu=Me/Mi;                                           % Ratio of Electron to ion mass used in coupling scheme
    gammai=3;                                          	% politropic constant ions(Isothermal csae only)


% presheath parameters:
    Jps=0;                                              % Total current density at preseath side
    ev=11604.9;
    K=1e3;
    Tpse=11*ev;                                  % Electron temperature at preseath side(k)
    Tpsi=10.5*ev;                                      % Ion temperature at preseath side(k)
    nps=1e20;                                           % Plasma density at preseath side(particle/m^3)
    cpsi=sqrt(kB*(Tpse+gammai*Tpsi)/Mi);                % Ion acoustic velocity at presheath side (m/s)
    upsi=-cpsi*cosd(Angle);                                         % Average fluid velocity along x-axis(sheath edge singularity condition)(m/s)

% Normalized quantities:
    JpsN=Jps/e/nps*sqrt(Mi/2/kB/Tpse);                  % Normalized current density
    TpsiN=Tpsi/Tpse;                                	% Normalized ion temperature
    upsiN=upsi/sqrt(2*kB*Tpse/Mi);                      % Normalized average fluid velocity

% sheath parameters:
    JL=Jps;                                             % Total current density at x=L(shaeth entrance)
    phiL=0;                                             % Potential at x=L(sheath entrance)
    rhoL=0;                                             % Total charge density at x=L(sheath entrance)
    L ='10*DLe';                                        % System length(length between the wall & sheath entarace)

% numerical input parameters:
    ntra=200;                                           % Total number of trajectories at sheath entarce between vicL & vimaxL
    phi0x='phi0*(1-X/L)';                               % Potential at X-position iss obtained by interpolation
    nx=50;                                              % Total number  of x-grids
    dt=1e-11;                                           % Time step size(s)
    deltaphi=1e-7;                                      % Desired accuracy in the iteration((phi2-phi1)>deltaphi)(volt)
    w=.03;                                              % Relaxitation parameter(relating to relaxitation scheme)
    nv=300;                                             % Total number of velocity grids

   phifNi=linspace(-5,-0.001,1001);                    % Initial guess (always negative)
   Tau_min=-5.8877;                              % Lower limit of variable 'Tau'. useless to keep below: -5.8877254338, since f -> 0.

   TcLiNi=linspace(-5.88,-0.001,90);                        % initial guess	%upper limit for TpsiN=.1 -> 2.2373356648402
   asymTcLiN=-4.5;                                     % limit below which we use the asymptotic expression
