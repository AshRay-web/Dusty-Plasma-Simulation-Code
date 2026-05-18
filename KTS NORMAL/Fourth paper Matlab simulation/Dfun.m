%specify the injection ion velocity distribution function for specified velocity vx
function df = Dfun(vx)                                     % Dfun is the name give to this file.(Defining the M.file as a function ,and is saved DFUN.m
global niL vmLi vtfi  vmaxLi vcLi                        % Importing global parameter from the mail file 
df=2*niL*exp(-(((vx-vmLi)^2)+(-5.6160e+005)^2+(-5.6160e+005)^2)/(vtfi)^2)/(vtfi*sqrt(pi))^3/(erf((-vmaxLi+vmLi)/vtfi)+erf((-vmLi+vcLi)/vtfi)); %Ion distribution function at x=L