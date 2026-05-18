     close all,                                              % Close all active figure,if any
     clear all,                                              % Clear all previous data,if any
     global niL vmLi vtfi kB vmaxLi vcLi Ef dx X nx          % Defiming global parameters


Input                                                        %% input('Type the name of your input file:'); 


NecessaryParameters                                          % Gives value of constant  parameters for sheath simulation by presheath-sheath coupling


     DLe=sqrt(epsilon*kB*Tfe/neL_/e^2);                      % Debye length of injected electron
     L=eval(L);                                              % Evaluate the value of L from
     vmaxLi=Tau_min*vtfi-vmLi;                               % Maximum  injection  velocity of ion  at x=L

I=1;                                                         % Initializing the number of iteration counter
     X=linspace(0,L,nx);                                     % Linear position grid with nx grid point between 0&L
     dx=X(2)-X(1);                                           % x-grid mesh size
     VLi=linspace(vmaxLi,vcLi,ntra);                         % Linear velocity grid with ntra point between vmaxLi & vcLi
     dVLi=abs(VLi(2)-VLi(1));                                % Velocity mesh grid size
if   dx/DLe > .8                                             % If Grid width is larger than Debye length ,stop the calculation 
     disp(' Attention: Your grid width is comparable/larger than Debye length.')
     disp('Press ENTER to continue or CTRL-C to stop.')
     pause
end                                                          %  Done the calcualtion if grid size is less then debye length
     phi=eval(phi0x); phi(1)=phi0; phi(nx)=phiL;             % Initial guess of potential profile(given at wall & sheath entrance) 
     
for  j=2:nx-1                                                % Selecting an inner jth x-grid point from 2 to nx-1.
     Ef(j)=.5*(phi(j-1)-phi(j+1))/dx;                        %  Calculate electric field at any grid point(j) in the simulation region
end                                                          % End of the loop started for j=2:nx-1
     Ef(1)=.5*(3*phi(1)-4*phi(2)+phi(3))/dx;                 % Electric field at(wall) x=0
     Ef(nx)=.5*(-3*phi(nx)+4*phi(nx-1)-phi(nx-2))/dx;        % Electric field at(sheath entrance) x=L
 DFAiondensity                                               % Calculate the ion density using DFA
     ne=neL_*exp(e*phi/kB/Tfe).*(1+erf(sqrt(e*(phi-phi(1))/kB/Tfe)));     % Electron density distribution at grid point(xj)
     ni=niL*ni/max(ni); ne=niL*ne/max(ne);                     % Check the density goes out of the boundary region
     nin=ni/niL;                                               % Normalized density of ion
     rho=Ze*ni-e*ne;                                           % Total charge density
     d=-2*ones(1,nx); d(1)=1; d(nx)=1; u=ones(1,nx-1);         % Diagonal value of the Matrix of Left hand side matrix
     ML=diag(d)+diag(u,1)+diag(u,-1); ML(1,2)=0;ML(nx,nx-1)=0; % Determinant of Mtrix of left hand side matrix
     clear d u                                                 % Clear the value of d & u
     poissionsolver                                            % Solve the poission equation and gives the new potential as 'phin'
     %plot(X/DLe,phi)
     %hold, pause(.1)
while max(abs(phin-phi))>deltaphi                              % Check the fluctuation of old & new potential
 I=I+1                                                         % Next iteration
      fluctuation=max(abs(phin-phi))                           % Shows the fluctuation of potential
      plot(X/DLe,phi)                                      % plot the graph between normalized distance and potential
      pause(.1)
            
      phi=phin;                                                % Replacing the old potential by calculated  new potential
  for j=2:nx-1                                                 % Selecting an inner jth x-grid point
      Ef(j)=.5*(phi(j-1)-phi(j+1))/dx;                         % Electric field at any j-th point in the simulation region
  end                                                          % End of loop started as for j=nx-2
      Ef(1)=.5*(3*phi(1)-4*phi(2)+phi(3))/dx;                  % Electric field at the x=0
      Ef(nx)=.5*(-3*phi(nx)+4*phi(nx-1)-phi(nx-2))/dx;         % Electric field at the x=L
      
  DFAiondensity                                                % Calculate the new ion density using DFA
  
      ne=neL_*exp(e*phi/kB/Tfe).*(1+erf(sqrt(e*(phi-phi(1))/kB/Tfe))); % Electron density distribution at grid point
      ni=niL*ni/max(ni); ne=niL*ne/max(ne);                    % Check the density goes out of the boundary region
      nin=ni/max(ni);                                          % Normalised ion density
      nen=ne/max(ne);                                          % Normalised electron density
      vcLe=sqrt(-2*e*phi0/Me);                                 % Cut off velocity of electron at x=L 
      ELe=0.5*Me*vtfe^2*(1.5-vcLe*De/vtfe/Ce/sqrt(pi))/e;
      Exe=0.5*Me*vtfe^2*(1.5-sqrt(e*(phi-phi(1))/kB/Tfe).*exp(-1*(e*(phi-phi(1))/kB/Tfe))./(1+erf(sqrt(e*(phi-phi(1))/kB/Tfe)))/sqrt(pi));
      Ee=Exe/e;
      Een=Exe/max(Exe);
      rho=Ze*ni-e*ne;                                          % Charge density distribution for new potential profile 
      poissionsolver                                           % Solving the poission equation for new potential
end                                                            % End the calculation while max(abs(phin-phi))>deltaphi
      plot(X/DLe,phi,'r')                                      % plot the graph between normalized distance and potential
      ne=neL_*exp(e*phi/kB/Tfe).*(1+erf(sqrt(e*(phi-phi(1))/kB/Tfe))); % Electron density distribution at grid point
      ni=niL*ni/max(ni); ne=niL*ne/max(ne);                    % Check the density goes out of the boundary region
      nin=ni/max(ni);                                          % Normalised ion density
      nen=ne/max(ne);       % Normalised electron density
      wi=ni*Vxm