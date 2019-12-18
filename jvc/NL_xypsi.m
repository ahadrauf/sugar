%function NL_angles2Fq(phi0,thetaB,psi0) 

function [x,y,psi,s,dx,i]=NL_xypsi(p,k,phib,psi0)
%By: jvclark - spring 2002

i=0;
resolution=100;
length=50e-6;

for theta=pi/2 : (((pi/2)-psi0)-(pi/2))/resolution : (pi/2)-psi0   
    i=i+1; %iterant
    
    phi=asin(sin(theta/2)/p); %
    h=2*p/k;
    x(i)=h*(cos(phib)-cos(phi));
    dx(i)=x(i) - length;

    %elliptic integration
    Ephib=quadl(@ellipticintegrand,0,phib,[],[],p,2); %elliptic integral of the second kind 
    Fphib=quadl(@ellipticintegrand,0,phib,[],[],p,1); %elliptic integral of the fist kind 
    Ephi=quadl(@ellipticintegrand,0,phi,[],[],p,2); %elliptic integral of the second kind 
    Fphi=quadl(@ellipticintegrand,0,phi,[],[],p,1); %elliptic integral of the first kind 
    
    y(i)=-(2*Ephib-Fphib-2*Ephi+Fphi)/k; %beam transverse deflection 
    s(i)=(Fphi-Fphib)/k; %s 
    psi(i)=pi/2 - 2*asin(p*sin(phi)); %angle along beam
end

function F=ellipticintegrand(x,p,kind)
%%By: jvclark - spring 2002
    if kind==1 %elliptic integral of the first kind
        F=1./sqrt(1-(p.*sin(x)).^2); %integrand of the first kind
    else %elliptic integral of the second kind
        F=sqrt(1-(p.*sin(x)).^2); %integrand of the second kind
    end


%[X,Y,RZ,S,DX]=NL_xypsi(p,k,phib,rz);figure(2);plot(X,Y,'g',S,Y,'b');x=X(11),y=Y(11),rz=RZ(11),dx=DX(11)

%PHI0=pi/2; THETAB=pi/2; PSI0=-pi/2.01; [x,y,rz,fx,fy,m,lm,p,k,phib]=NL_angles2Fq2(pi/2,pi/2,PSI0);dx=x,y,rz,  [X,Y,RZ,S,DX]=NL_xypsi(p,k,phib,rz);figure(2);plot(X,Y,'g',S,Y,'b');dx=DX(11),y=Y(11),rz=RZ(11), title('beam curves');xlabel('green=>y(x),  blue=>y(s)');

%PHI0=pi/2; THETAB=pi/2; PSI0=-pi/2.0000001; [x,y,rz,fx,fy,m,lm,p,k,phib]=NL_angles2Fq2(pi/2,pi/2,PSI0);dx=x,y,rz,  [X,Y,RZ,S,DX,i]=NL_xypsi(p,k,phib,rz);figure(2);plot(X,Y,'g',S,Y,'b');dx=DX(i),y=Y(i),rz=RZ(i), title('Beam Curves. green=y(x), blue=y(s), L=50um, h=w=2um, E=165e9Pa');string=strcat('End values: dx=',real(num2str(dx,3)),'m, y=',real(num2str(y,3)),'m, rz=',real(num2str(rz,3)),'rad, s=',ral(num2str(S(i),3)),'m');xlabel(string);string=strcat('Fx=',real(num2str(fx,3)),'N, Fy=',real(num2str(fy,3)),'N, Mz=',real(num2str(m,3)),'Nm');ylabel(string);
%PHI0=pi/2; THETAB=pi/2; PSI0=pi/1.6; [x,y,rz,fx,fy,m,lm,p,k,phib]=NL_angles2Fq2(pi/2,pi/2,PSI0);dx=x,y,rz,  [X,Y,RZ,S,DX,i]=NL_xypsi(p,k,phib,rz);figure(2);plot(X,Y,'g',S,Y,'b');dx=DX(i),y=Y(i),rz=RZ(i), title('Beam Curves. green=y(x), blue=y(s), L=50um, h=w=2um, E=165e9Pa');string=strcat('End values: dx=',real(num2str(dx,3)),'m, y=',real(num2str(y,3)),'m, rz=',real(num2str(rz,3)),'rad, s=',real(num2str(S(i),3)),'m');xlabel(string);string=strcat('Fx=',real(num2str(fx,3)),'N, Fy=',real(num2str(fy,3)),'N, Mz=',real(num2str(m,3)),'Nm');ylabel(string);
