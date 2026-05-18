
%Input file that is necessary to run the simulation.
%2025-11-26
%Input parameters in SI units
angle=30;
theta=(pi/180)*angle;
B=0.05;
kB=1.3806e-23;
e=1.6022e-19;
epsilon_not=8.8542e-12;
me=9.1094e-31;
mi=1.6737e-27;
mu=me/mi;
gammai=5/3;
gammae=1;


% constants
eV   = 1.16045221e4;        % eV in Kelvin equivalent
Jps=0;
Tpse = 10*eV;               % electron temperature (K)
Tpsi = 1*eV;          % ion temperature (K)
nps  = 1e20;                % reference plasma density [m^-3]

% --- plasma densities ---
niL  = nps; 

% --- dust parameters ---
Zd   = 1e3;                 % dust charge number
Qd0  = -Zd;                 % normalized dust charge
nd = 1e14;
neL   = niL+nd*Qd0 ;


% --- ion acoustic speed ---
cpsi = sqrt((kB*Tpse + gammai*kB*Tpsi)/mi);
upsi = -cpsi;


JpsN=Jps/e/nps/sqrt(2*kB*Tpse/mi);
TpsiN=Tpsi/Tpse;
upsiN=upsi/sqrt(2*kB*Tpse/mi);

JL=Jps;
phiL=0;
rhoL=0;
L='10*DLe';

ntra=200;                                   % Total number of trajectories at sheath entrance between vicL & vimaxL
phiwx='phiw*(1-X/L)';                       % Potential at X-position is obtained by interpolation
nx=50;                                     % Total number of x-grids

dt=0.5e-12;                                    % Time step size(s)
deltaphi=1e-7;                              % Desired accuracy in the iteration ((phi2-phi1)>deltaphi)(volt)
w=0.08;                                     % Relaxitation parameter(relating to relaxitation scheme)
nv=300;                                     % Total number of velocity grids
phiNi=linspace(-5,-0.001,1001);             % Initial guess (always negative)
Tau_min=-6;
TcLiNi=linspace(-6,0,1001);
asymTcLiN=-4.5;

