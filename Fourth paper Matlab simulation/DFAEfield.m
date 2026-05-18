%prepared  at 05/10/2012 by Roshan Chalise,cdp
function field=DFAEfield(x)
global Ef dx X nx
if x==0
   field=Ef(1);
elseif x==X(nx)
   field=Ef(nx);
else
   j=ceil(x/dx);
   field=(Ef(j)*(X(j+1)-x)+Ef(j+1)*(x-X(j)))/dx;
end  