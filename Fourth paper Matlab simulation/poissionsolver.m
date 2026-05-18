% Prepared at 05/09/2012 by Roshan Chalise,CDP
% solves the poissons rquation for given charge density and potential at
% the boundaries fixed to some constant value
MR(1)=phi(1); MR(nx)=phi(nx); MR(2:nx-1)=-dx^2*rho(2:nx-1)/epsilon;   %RHS Matrix of poission solution
phin=ML\MR'; phin=phin';                                             % Solution of 'ML*phin=MR' Ml is defined in main command file
phin= w*phin + (1-w)*phi;                                            % Relaxation scheme
phin(1)=phi(1);                                                       % Potential at the wall(x=0)
phin(nx)=phi(nx);                                                     % Potential at injection point (x=L)
