%Distribution function of ions
function df = Dfun(vx,vy,vz)
global niL vmLi vti vcLi
df =(2*niL/(vti*sqrt(pi))^3)/(1+erf((vcLi-vmLi)/vti))*exp(-(((vx-vmLi)^2+(vy)^2+(vz)^2)/(vti^2)));
