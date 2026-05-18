% Coupling sheath with presheath.
% For given presheath plasma parameters (nps, Tpsi, Tpse, Jps, gammai, mu, etc.) this program
% calculates sheath plasma parameters (Ai, Ae, vcLi, vmLi, Tfe, Tfi, phi0, etc.) required for
% our 1d1v plasma-sheath simulation.
%
% This proceeds as following:
%
% 1) Provide the name of your input file, where all required input parameters are specified.
%
% 2) Solves the electron irreducible equation to obtain 'phifN', and then the other corresponding
% electron parameters (TfeN, AeN, phi0N, etc.) are calculated.
%
% 3) Solves the ion irreducible equation to obtain 'TaucLiN', and then the other corresponding0
% ion parameters (TfiN, AiN, vmLiN, vcLiN, etc.) are calculated.
%
% 4) Once all normalized parameters are known, the dimensionless physical parameters are calculated
%    using their respective normalizing equations.
%



% solving the electron irreducible equation for 'phifN' --- start --- :
%
% Calculating the lhs (Le) of the equation (independent of 'JpsN' & 'TpsiN'):
for rphifN=1:length(phifNi)
   phifN=phifNi(rphifN);
	Ce=1+erf(sqrt(-phifN));
	De=exp(phifN);
	TfeN=1/abs(1 - sqrt(-4*phifN/pi)*De/Ce - 2/pi*(De/Ce)^2);
   Le(rphifN)=De/Ce*sqrt(TfeN/pi);		% lhs of the electron irreducible equation
end	%for r=1:length(phifNi)
%      
% Calculating the rhs (Re) of the equation (function of 'JpsN', 'TpsiN', 'gammai' and 'mu'):
Re=sqrt(mu)*(JpsN+sqrt(.5*(1+gammai*TpsiN)));	% rhs of the equation
phifN_sol=interp1(Le, phifNi, Re);		% Interpolating for the solution (phifN) for which Le = Re.
%      
% Thus obtained 'phifN_sol' is the required value.
if isnan(phifN_sol)
   disp('Could not find the solution. Change the input-range for phifN and start again.')
   pause
else
   phifN=phifN_sol
	% Calculating other parameters:
	Ce=1+erf(sqrt(-phifN));
	De=exp(phifN);
	TfeN=1/abs(1 - sqrt(-4*phifN/pi)*De/Ce - 2/pi*(De/Ce)^2);
	AeN=1/Ce/sqrt(TfeN);
end		% if isnan(phifN_sol)
% solving the electron irreducible equation for 'phifN' --- end ---


% solving the ion irreducible equation for 'TcLiN' ------------- start ------------
%
% Calculating the rhs (Ri) of the equation (function of 'phifN'; independent of 'JpsN' & 'TpsiN'):
Ri=(sqrt(pi) + De/Ce/sqrt(-phifN))/TfeN;
%
% Calculating the lhs (Li) of the equation:
rr=1;
for r=1:length(TcLiNi)
   TcLiN=TcLiNi(r);
   if TcLiN<asymTcLiN
      % using asymptotic expansion: erf(x)=1-exp(-x^2)/sqrt(pi)/x*(1 - 1/2/x^2 + 3/2^2/x^4 - 1.3.5/2^3/x^6 + ...)
      Ci=exp(-TcLiN^2)/abs(TcLiN)/sqrt(pi)*(1 - .5/TcLiN^2 + .75/TcLiN^4 - 15/8/TcLiN^6 + 105/16/TcLiN^8 - ...
         945/32/TcLiN^10 + 10395/64/TcLiN^12 - 135135/128/TcLiN^14 + 2027025/256/TcLiN^16 - 34459425/512/TcLiN^18 + ...
         654729075/1024/TcLiN^20 - 13749310575/2048/TcLiN^22 + 316234143225/4096/TcLiN^24 - 7905853580625/8192/TcLiN^26);
   else	%if TcLiN<asymTcLiN
      Ci=1+erf(TcLiN);
   end	%if TcLiN<asymTcLiN
   Di=exp(-TcLiN^2);
   TfiN=TpsiN/abs(1 - 2*TcLiN*Di/Ci/sqrt(pi) - 2*(Di/Ci)^2/pi);
   vmLiN=-sqrt(0.5*(1+gammai*TpsiN)) + Di/Ci*sqrt(TfiN/pi);
   a=vmLiN/sqrt(TfiN);	% abbrebiated for simplicity
   if TcLiN < -a	% checking the integrability
      % Integratinng:
      Tau=linspace(Tau_min,TcLiN,ntra);	% Discretizing for integration
      dTau=Tau(2) - Tau(1);					% Width of the Tau-grid
      for rTau=1:ntra
         FnTau(rTau)=exp(-Tau(rTau)^2) / (Tau(rTau) + a)^2;
      end
      LHS = dTau * ( sum(FnTau) - 0.5*(FnTau(1) +FnTau(ntra)) );
      Li(rr) = LHS/Ci/TfiN;
      TcLiNM(rr)=TcLiN;				% storing the values for which the equation is integrable. 
      rr=rr+1;
   end	% if TcLiN < -a
end	%for r=1:length(phifNi)
%
% obtaining the solution by locating the point of intersection of 'Li' & 'Ri':
%%plot(TcLiNM, Li)
%%hold
%%plot([min(TcLiNM) max(TcLiNM)],[Ri Ri],'k--')
TcLiN_sol=interp1(Li,TcLiNM,Ri);	% required solution
%
% This method of obtaining the solution can be RISKY. The point where we have found 'Li=Ri' is
% the point of marginal validity of the Bohm's criterion, which, in our case, is Li <= Ri. When
% we consider the marginal point it may be possible that because of numerical limitations it
% lies in fact in the region, where the Bohm's criterion is not satisfied. Hence, in order to
% be sure we take a point very close to the zero-point but still lying in the negative side of
% sure side of (Li - Ri), so that we may say that Li <~ Ri, in place of Li <= Ri.
% We have tested it many times and concluded that it has also negligible effect to the other
% ion parameters derived therefrom.
a=find(Li-Ri<0);						% points where Li < Ri
a=a(length(a));						% points closest to the zero-point but still Li < Ri
TcLiN_sol=TcLiNM(a);					% value of the closest point to the zero-point => Li <~ Ri
%%plot(TcLiNM(a), Li(a),'mo')
% solving the ion irreducible equation for 'TcLiN' ------------- end ------------


% calculating other dimensionless ion parameters:
%
% Thus obtained 'TcLiN_sol' is the required value.
if isnan(TcLiN_sol)
   disp('Could not find the solution. Change the input-range for TcLiN and start again.')
   pause
else
   TcLiN=TcLiN_sol
	% Calculating other parameters:
   if TcLiN<asymTcLiN
      % using asymptotic expansion: erf(x)=1-exp(-x^2)/sqrt(pi)/x*(1 - 1/2/x^2 + 3/2^2/x^4 - 1.3.5/2^3/x^6 + ...)
      Ci=exp(-TcLiN^2)/abs(TcLiN)/sqrt(pi)*(1 - .5/TcLiN^2 + .75/TcLiN^4 - 15/8/TcLiN^6 + 105/16/TcLiN^8 - ...
         945/32/TcLiN^10 + 10395/64/TcLiN^12 - 135135/128/TcLiN^14 + 2027025/256/TcLiN^16 - 34459425/512/TcLiN^18 + ...
         654729075/1024/TcLiN^20 - 13749310575/2048/TcLiN^22 + 316234143225/4096/TcLiN^24 - 7905853580625/8192/TcLiN^26);
   else	%if TcLiN<asymTcLiN
      Ci=1+erf(TcLiN);
   end	%if TcLiN<asymTcLiN
   Di=exp(-TcLiN^2);
   TfiN=TpsiN/abs(1 - 2*TcLiN*Di/Ci/sqrt(pi) - 2*(Di/Ci)^2/pi);
   vmLiN=-sqrt(0.5*(1+gammai*TpsiN)) + Di/Ci*sqrt(TfiN/pi);
   AiN=1/Ci/sqrt(TfiN);
	vcLiN=vmLiN + TcLiN*sqrt(TfiN);
end	%if isnan(TcLiN_sol)

% calculating dimensional physical parameters:
%
% electron dimensional parameters:
Tfe=TfeN*Tpse;
vtfe=sqrt(2*kB*Tfe/Me);
Ae=AeN*nps*sqrt(2*Me/pi/kB/Tpse);
neL_=Ae*vtfe*sqrt(pi)/2;
phi0=phifN*kB*Tfe/e;
%
% ion dimensional parameters:
Tfi=TfiN*Tpse;
vtfi=sqrt(2*kB*Tfi/Mi);
Ai=AiN*nps*sqrt(2*Mi/pi/kB/Tpse);
niL=nps;
vmLi=vmLiN*sqrt(2*kB*Tpse/Mi);
vcLi=vcLiN*sqrt(2*kB*Tpse/Mi);


