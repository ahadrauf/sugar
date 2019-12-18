%function NL_angles2Fq(phi0,thetaB,psi0) 

function [x,y,rz,Fphi0-Fphib,Fx,Fy,M,L3,p,k,phib,rz]=NL_integrations1(phi0,thetaB,psi0) 
%By: jvclark - spring 2002
    
%parameters
    theta0=thetaB + psi0; %thetaB => angle of force, psi0 => tip angle rz.
    p=sin(theta0/2); %elliptic integral modulus
    L3k=2*p*cos(phi0); %associated with moment 
    phib=asin(sin(thetaB/2)/p); %angle associated with beam length and moment
    
%elliptic integration
    Ephib=quadl(@ellipticintegrand,0,phib,[],[],p,2); %elliptic integral of the second kind 
    Fphib=quadl(@ellipticintegrand,0,phib,[],[],p,1); %elliptic integral of the fist kind 
    Ephi0=quadl(@ellipticintegrand,0,phi0,[],[],p,2); %elliptic integral of the second kind 
    Fphi0=quadl(@ellipticintegrand,0,phi0,[],[],p,1); %elliptic integral of the first kind 

    L1minusXoverL1=2*p*(cos(phib)-cos(phi0))/L1k; %beam shortening 
    YoverL1=(2*Ephib-Fphib)/L1k - (2*Ephi0-Fphi0)/L1k; %beam transverse deflection 
    FL2EIx=L1k^2*cos(thetaB); %nondimensional force
    FL2EIy=L1k^2*sin(thetaB); %nondimensional force
    MLEI=L3k*L1k; %nondimensional moment
    psi0overpi2=psi0/(pi/2); %nondimensional tip angle
 
    %re-sign to Sugar directions
        %L1minusXoverL1=-1+L1minusXoverL1;
        %YoverL1=-YoverL1;
        %psi0overpi2=-psi0overpi2;
        %FL2EIx=-FL2EIx;
        %FL2EIy=-FL2EIy;
        %MLEI=-MLEI;
    
    x=L1minusXoverL1*L1-L1;
    y=YoverL1*-L1;
    Fx=-P*cos(thetaB);
    Fy=-P*sin(thetaB);
    M=L3k*Fy/k;
    rz=-psi0; %change back
    L3=L3k/k;
    
%clear all;tic;i=0;for frac=0.01:0.01:1;i=i+1;phi0(i)=pi/2 *frac;[X(i),Y(i),PSI(i),FX(i),FY(i),M(i),L1(i),L3(i),x(i),y(i),rz(i),fx(i),fy(i),m(i),lm(i)]=NL_angles2Fq(phi0(i),pi/2,-pi/2 *0.5);end;toc;figure(1);phi0=phi0/(pi/2);phi0=L3;plot(phi0,FY,'y',phi0,X,'g',phi0,Y,'b',phi0,PSI,'r',phi0,FX,'m',phi0,M,'c');xlabel('xaxis: phi0/(pi/2).   yaxis: green=X  blue=Y  red=Psi0  magenta=Fx  yellow=Fy  cyan=M');title('sweep moment arm => phi0 => L3');grid on;figure(2);plot(lm,x,'g',lm,y,'b',lm,rz,'r',lm,fx,'m',lm,fy,'y',lm,m,'c');xlabel('xaxis: phi0/(pi/2).   yaxis: green=X  blue=Y  red=Psi0  magenta=Fx  yellow=Fy  cyan=M');title('sweep moment arm => L3');grid on;
%clear all;tic;i=0;for frac=0.01:0.01:1;i=i+1;phi0(i)=pi/2 *frac;[X(i),Y(i),PSI(i),FX(i),FY(i),M(i),L1(i),L3(i),x(i),y(i),rz(i),fx(i),fy(i),m(i),lm(i)]=NL_angles2Fq(phi0(i),pi/2,-pi/2 *0.5);end;toc;figure(1);phi0=phi0/(pi/2);phi0=L3;plot(phi0,FY,'y',phi0,X,'g',phi0,Y,'b',phi0,PSI,'r',phi0,FX,'m',phi0,M,'c');xlabel('xaxis: phi0/(pi/2).   yaxis: green=X  blue=Y  red=Psi0  magenta=Fx  yellow=Fy  cyan=M');title('sweep moment arm => phi0 => L3');grid on;figure(2);plot(lm,x,'g',lm,y,'b',lm,fy,'y',lm,m,'c');xlabel('xaxis: phi0/(pi/2).   yaxis: green=X  blue=Y  red=Psi0  magenta=Fx  yellow=Fy  cyan=M');title('sweep moment arm => L3');grid on;
        
function F=ellipticintegrand(x,p,kind)
%%By: jvclark - spring 2002
    if kind==1 %elliptic integral of the first kind
        F=1./sqrt(1-(p.*sin(x)).^2); %integrand of the first kind
    else %elliptic integral of the second kind
        F=sqrt(1-(p.*sin(x)).^2); %integrand of the second kind
    end


%clear all;tic;i=0;for frac=0.01:0.01:1;i=i+1;phi0(i)=pi/2 *frac;[x(i),y(i),rz(i),fx(i),fy(i),m(i),lm(i)]=NL_angles2Fq2(phi0(i),pi/2,-pi/2 *0.5);end;toc;figure(1);plot(lm,x,'g',lm,y,'b',lm,fy,'y',lm,m,'c');xlabel('xaxis: phi0/(pi/2).   yaxis: green=X  blue=Y  red=Psi0  magenta=Fx  yellow=Fy  cyan=M');title('sweep moment arm => L3');grid on;  
%clear all;tic;i=0;for phi=pi/2:-0.01:pi/8; i=i+1; phi0(i)=phi;[x(i),y(i),rz(i),fx(i),fy(i),m(i),lm(i)]=NL_angles2Fq2(phi0(i),pi/2,-pi/2 *0.5);end;toc;figure(1);X=phi0;plot(X,x,'g',X,y,'b',X,fy,'y',X,m*1e5,'c');xlabel('xaxis: phi0/(pi/2).   yaxis: green=X  blue=Y  red=Psi0  magenta=Fx  yellow=Fy  cyan=M');title('sweep moment arm => L3');grid on;

%sweep psi0
%clear all;tic;i=0;for psi=pi/2*0.9:-0.01:-pi/2*0.9; i=i+1; psi0(i)=psi;[x(i),y(i),rz(i),fx(i),fy(i),m(i),lm(i)]=NL_angles2Fq2(pi/2,pi/2,psi0(i));end;toc;figure(1);X=psi0;plot(X,x,'g',X,y,'b',X,fy,'y',X,m*1e4,'c',X,fx,'m');xlabel('xaxis: psi0.   yaxis: green=X  blue=Y  red=Psi0  magenta=Fx  yellow=Fy  cyan=M');title('sweep psi0=pi/2 to -pi/2,    phi0=pi/2,    thetaB=pi/2');grid on;
% clear all;tic;i=0;for psi=pi/2*0.9:-0.01:-pi/2*0.9; i=i+1; psi0(i)=psi;[x(i),y(i),rz(i),fx(i),fy(i),m(i),lm(i),p(i),k(i),phib(i)]=NL_angles2Fq2(pi/2,pi/4,psi0(i));end;toc;figure(1);X=psi0/(pi/2);plot(X,x,'g',X,y,'b',X,fy,'y',X,m*1e4,'c',X,fx,'m');xlabel('xaxis: psi0/(pi/2).   yaxis: green=X  blue=Y  red=Psi0  magenta=Fx  yellow=Fy  cyan=M');title('sweep psi0=pi/2 to -pi/2,    phi0=pi/2,    thetaB=pi/4');grid on;
