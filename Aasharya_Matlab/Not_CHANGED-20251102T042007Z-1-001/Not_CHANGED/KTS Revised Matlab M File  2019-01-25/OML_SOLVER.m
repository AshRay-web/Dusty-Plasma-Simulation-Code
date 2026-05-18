%close all;
%clear all;

vi=sqrt(vxi(200,:).^2+vyi(200,:).^2+vzi(200,:).^2);
rd=10e-6;
Cs=sqrt(kB*Te/mi);                             
b=1;
a=sqrt(8*kB*Te/pi/me)/Cs;
%==========================================

Qdn0=-Zd;
Qd_r=zeros(nx,1);
%===================================================
for k=1:nx
    Ne=ne(k); % Not divided by nps as we already provide the normalized one
    Ni=ni(k);
    ui=abs(vi(k))/Cs; 
%===========For Qd less than zero.
 Fqdn = @(Qdn) Ne*a*exp(Qdn)-Ni*ui*(1-(2*b*Qdn)/ui^2);
 dFqdn = @(Qdn) Ne*a*exp(Qdn)+2*b*Ni/ui;
Qdn1=Qdn0-Fqdn(Qdn0)/dFqdn(Qdn0);     %%%%%
error_n=abs(Qdn1-Qdn0);
       while  error_n > epsilon_not
        Qdn1=Qdn0-Fqdn(Qdn0)/dFqdn(Qdn0);
        error_n=abs(Qdn1-Qdn0);
        Qdn0=Qdn1;
       end
  Qdn(k)=Qdn0;
end 
C0=4*pi*epsilon_not;

qdn=C0*rd*kB*Te*Qdn/e;

phid=qdn/C0/rd;
