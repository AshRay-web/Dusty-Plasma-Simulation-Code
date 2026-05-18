MR(1)=phi(1);
MR(nx)=phi(nx);
MR(2:nx-1)=-dx^2*rho(2:nx-1)/epsilon_not;
phin=ML\MR';
phin=phin';
phin=w*phin+(1-w)*phi;
phin(1)=phi(1);
phin(nx)=phi(nx);