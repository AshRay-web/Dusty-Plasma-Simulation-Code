%This file is used to calculate the ion density.
Q=e*dt/mi;
for T=1:ntra
    m=1;
    vxL=vxLi(T);
    vyL=vyLi(T);
    vzL=vzLi(T);
    vxm=vxL;
    vym=vyL;
    vzm=vzL;
    vxmplusonefourth_time=vxm+0.25*Q*Efieldsolver(L)-0.25*Q*B*vzm*sin(theta);
    Xmplushalf=L+0.5*dt*vxmplusonefourth_time;
    vymplusonefourth_time=vym+0.25*Q*B*vzm*cos(theta);
    vzmplushalf=vzm+0.5*Q*B*(vxmplusonefourth_time*sin(theta)-vymplusonefourth_time*cos(theta));
    
    vxm(m)=vxm;
    vym(m)=vym;
    vzm(m)=vzm;
    Xm2(m)=Xmplushalf;
    vxm2(m)=vxm;
    vym2(m)=vym;
    vzm2(m)=vzm;
    if Xm2(m)<=X(nx-1);
        disp('time step is too large');
        pause
    end
    while Xm2(m)>0 & Xm2(m)<=L
        for m=m+1
            vxm(m)=vxm(m-1)+Q*Efieldsolver(Xm2(m-1))-Q*B*vzm(m-1)*sin(theta);
            vym(m)=vym(m-1)-Q*B*vzm(m-1)*cos(theta);
            vxmminushalf=0.5*(vxm(m)+vxm(m-1));
            vymminushalf=0.5*(vym(m)+vym(m-1));
            vzm(m)=vzm(m-1)+Q*B*(vxm(m-1)*sin(theta)-vym(m-1)*cos(theta));
            Xm2(m)=Xm2(m-1)+dt*vxm(m);
        end
    end
            for r=1:m-1
            vxm2(r)=0.5*(vxm(r)+vxm(r+1));
            vym2(r)=0.5*(vym(r)+vym(r+1));
            vzm2(r)=0.5*(vzm(r)+vzm(r+1));
    end
    Xt=[L,Xm2(1:m-1)];
    vxt=[vxL,vxm2];
    vyt=[vyL,vym2];
    vzt=[vzL,vzm2];
    a=((Xt(m-2)-Xt(m))^2*(vxt(m-1)-vxt(m-2))-(Xt(m-2)-Xt(m-1))^2*(vxt(m)-vxt(m-2)))/(Xt(m-2)-Xt(m-1))/(Xt(m-2)-Xt(m))/(Xt(m-1)-Xt(m));
    b=((Xt(m-2)-Xt(m))*(vxt(m-1)-vxt(m-2))-(Xt(m-2)-Xt(m-1))*(vxt(m)-vxt(m-2)))/(Xt(m-2)-Xt(m))/(Xt(m-2)-Xt(m-1))/(Xt(m)-Xt(m-1));
    vx0=vxt(m-2)+a*Xt(m-2)+b*Xt(m-2)^2;
    
    a=((Xt(m-2)-Xt(m))^2*(vyt(m-1)-vyt(m-2))-(Xt(m-2)-Xt(m-1))^2*(vyt(m)-vyt(m-2)))/(Xt(m-2)-Xt(m-1))/(Xt(m-2)-Xt(m))/(Xt(m-1)-Xt(m));
    b=((Xt(m-2)-Xt(m))*(vyt(m-1)-vyt(m-2))-(Xt(m-2)-Xt(m-1))*(vyt(m)-vyt(m-2)))/(Xt(m-2)-Xt(m))/(Xt(m-2)-Xt(m-1))/(Xt(m)-Xt(m-1));
    vy0=vyt(m-2)+a*Xt(m-2)+b*Xt(m-2)^2;
    
    a=((Xt(m-2)-Xt(m))^2*(vzt(m-1)-vzt(m-2))-(Xt(m-2)-Xt(m-1))^2*(vzt(m)-vzt(m-2)))/(Xt(m-2)-Xt(m-1))/(Xt(m-2)-Xt(m))/(Xt(m-1)-Xt(m));
    b=((Xt(m-2)-Xt(m))*(vzt(m-1)-vzt(m-2))-(Xt(m-2)-Xt(m-1))*(vzt(m)-vzt(m-2)))/(Xt(m-2)-Xt(m))/(Xt(m-2)-Xt(m-1))/(Xt(m)-Xt(m-1));
    vz0=vzt(m-2)+a*Xt(m-2)+b*Xt(m-2)^2;
    
    Xt=[Xt,0]; 
    vxt=[vxt,vx0];
    vyt=[vyt,vy0];
    vzt=[vzt,vz0];
    vxti=interp1(Xt,vxt,X);
    vyti=interp1(Xt,vyt,X);
    vzti=interp1(Xt,vzt,X);
    vxi(T,:)=vxti;
    vyi(T,:)=vyti;
    vzi(T,:)=vzti;
    Df(T)=Dfun(vxL,vyL,vzL);
    clear vxm vym vzm vxm2 vym2 vzm2 vxt vyt vzt vxti vyti vzti Xm2 Xt a b r m vx0 vyo vzo
end 
ni=zeros(1,nx);
for j=1:nx
    v_xi=vxi(:,j); 
    v_yi=vyi(:,j); 
    v_zi=vzi(:,j);
    for T=1:ntra-1;
        ni(j)=ni(j)+0.5*(Df(T)+Df(T+1))*abs(v_xi(T+1)-v_xi(T));
    end
end