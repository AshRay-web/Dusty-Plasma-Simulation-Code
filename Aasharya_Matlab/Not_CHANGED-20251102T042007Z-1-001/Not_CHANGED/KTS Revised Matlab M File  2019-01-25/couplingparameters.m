
%necessary parameters for simulation
for rphiN=1:length(phiNi);
    phiN=phiNi(rphiN);
    Ce=1+erf(sqrt(-phiN));
    De=exp(phiN);
    TeN=1/abs(1-sqrt(-4*phiN/9/pi)*De/Ce-(2/3/pi)*(De/Ce)^2);
    Le(rphiN)=(De/Ce)*sqrt(TeN/pi);
end
Re=sqrt(mu)*(JpsN+sqrt(0.5*(1+gammai*TpsiN)));
phiN_sol=interp1(Le,phiNi, Re);
if isnan(phiN_sol)
    disp('could not find solution')
    pause
else
    phiN=phiN_sol;
    Ce=1+erf(sqrt(-phiN));
    De=exp(phiN);
    TeN=1/abs(1-sqrt(-4*phiN/9/pi)*De/Ce-(2/3/pi)*(De/Ce)^2);
    AeN=1/Ce/4/(TeN)^1.5;
end
%Ri=(sqrt(pi)+De/Ce/sqrt(-phiN))/TeN;

Ri = ((sqrt(pi) + De/Ce/sqrt(-phiN)) / TeN)*(1 + (nd/niL) * Qd0 ); %Modifiled
disp(Qd0)
rr=1;
for r=1:length(TcLiNi);
    TcLiN=TcLiNi(r);
    if TcLiN<asymTcLiN
        Ci=2-exp(-TcLiN^2)/(TcLiN)/sqrt(pi)*(1-0.5/TcLiN^2+3/(2*TcLiN^2)^2-15/(2*TcLiN^2)^3+105/(2*TcLiN^2)^4-945/(2*TcLiN^2)^5+10395/(2*TcLiN^2)^6-135135/(2*TcLiN^2)^7+2027025/(2*TcLiN^2)^8-34459425/(2*TcLiN^2)^9+654729075/(2*TcLiN^2)^10-13749310575/(2*TcLiN^2)^11+316234143225/(2*TcLiN^2)^12-7905853580625/(2*TcLiN^2)^13);
        %Ci=exp(-TcLiN^2)/abs(TcLiN)/sqrt(pi)*(1-0.5/TcLiN^2+3/(2*TcLiN^2)^2-15/(2*TcLiN^2)^3+105/(2*TcLiN^2)^4-945/(2*TcLiN^2)^5+10395/(2*TcLiN^2)^6-135135/(2*TcLiN^2)^7+2027025/(2*TcLiN^2)^8-34459425/(2*TcLiN^2)^9+654729075/(2*TcLiN^2)^10-13749310575/(2*TcLiN^2)^11+316234143225/(2*TcLiN^2)^12-7905853580625/(2*TcLiN^2)^13);
    else
        Ci=1+erf(TcLiN);
    end
    Di=exp(-TcLiN^2);
    TiN=TpsiN/abs(1-2*TcLiN*Di/Ci/sqrt(9*pi)-(2/3/pi)*(Di/Ci)^2);
    vmLiN=-sqrt(0.5*(1+gammai*TpsiN))+Di/Ci*sqrt(TiN/pi);
    a=vmLiN/sqrt(TiN);
    if TcLiN<-a;
        Tau=linspace(Tau_min, TcLiN,ntra);
        dTau=Tau(2)-Tau(1);
        for rTau=1:ntra;
            FnTau(rTau)=exp(-Tau(rTau)^2)/(Tau(rTau)+a)^2;
        end
        LHS=dTau*(sum(FnTau)-0.5*(FnTau(1)+FnTau(ntra)));
        Li(rr)=LHS/Ci/TiN;
        TcLiNM(rr)=TcLiN;
        rr=rr+1;
    end 
end
TcLiN_sol=interp1(Li,TcLiNM,Ri);
a=find(Li-Ri<0);
a=a(length(a));
TcLiN_sol=TcLiNM(a);
if isnan(TcLiN_sol);
    disp('could not find solution');
    pause
else 
    TcLiN=TcLiN_sol;
    if TcLiN<asymTcLiN;
        Ci=2-exp(-TcLiN^2)/(TcLiN)/sqrt(pi)*(1-0.5/TcLiN^2+3/(2*TcLiN^2)^2-15/(2*TcLiN^2)^3+105/(2*TcLiN^2)^4-945/(2*TcLiN^2)^5+10395/(2*TcLiN^2)^6-135135/(2*TcLiN^2)^7+2027025/(2*TcLiN^2)^8-34459425/(2*TcLiN^2)^9+654729075/(2*TcLiN^2)^10-13749310575/(2*TcLiN^2)^11+316234143225/(2*TcLiN^2)^12-7905853580625/(2*TcLiN^2)^13);
        %Ci=exp(-TcLiN^2)/abs(TcLiN)/sqrt(pi)*(1-0.5/TcLiN^2+3/(2*TcLiN^2)^2-15/(2*TcLiN^2)^3+105/(2*TcLiN^2)^4-945/(2*TcLiN^2)^5+10395/(2*TcLiN^2)^6-135135/(2*TcLiN^2)^7+2027025/(2*TcLiN^2)^8-34459425/(2*TcLiN^2)^9+654729075/(2*TcLiN^2)^10-13749310575/(2*TcLiN^2)^11+316234143225/(2*TcLiN^2)^12-7905853580625/(2*TcLiN^2)^13);
    else
     Ci=1+erf(TcLiN);
    end
    Di=exp(-TcLiN^2);
    TiN=TpsiN/abs(1-2*TcLiN*Di/Ci/sqrt(9*pi)-(2/3/pi)*(Di/Ci)^2);
    vmLiN=-sqrt(0.5*(1+gammai*TpsiN))+Di/Ci*sqrt(TiN/pi);
    AiN=1/Ci/4/(TiN)^1.5;
    vcLiN=vmLiN+TcLiN*sqrt(TiN);
end
Te=TeN*Tpse;
vte=sqrt(2*kB*Te/me);

niL  = nps;                 % assume ion density
             
neL   = niL+nd*Qd0 ;        % electron density

Ae=AeN*nps/(pi*kB*Tpse/2/me)^1.5;
phiw=phiN*kB*Te/e;
Cs = sqrt(kB*Te/mi);

%for ion sheath parameters
Ti=TiN*Tpse;
vti=sqrt(2*kB*Ti/mi);
Ai=AiN*nps/(pi*kB*Tpse/2/mi)^1.5;
vmLi=vmLiN*sqrt(2*kB*Tpse/mi);
vcLi=vcLiN*sqrt(2*kB*Tpse/mi);
