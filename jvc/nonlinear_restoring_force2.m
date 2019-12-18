%function Z=nonlineardeflection(x,y,rz)
%Given x, y, and rz displacements, 
%determine angles thetaB and phi0 (see derivation for term identification).

function [phi0,thetaB,psi0,x,y,rz,fx,fy,m,F,closeness,p_,k_,phib_]=nonlinear_restoring_force2(x_orig,y_orig,rz_orig) 
%jvclark - May 2002 

%global vars ONLY for the functions below
global L1k L3k L1minusXoverL1 YoverL1 closeness k P E I FL2EIx FL2EIy MLEI p phib 

%sign switch for formulation
    psi=-rz_orig; 
    x=-x_orig;
    y=-y_orig;

%material     
    E=165e9;
    I=(2e-6)^4/12;
    L1=50e-6;

%starting guess
    phi0 = pi/2 *0.8;
    thetaB = pi/2 *0.6;
    psi0 = psi; %external input
    
%true angles from minimization
    v=fminsearch(@nonlineardeflection,[phi0,thetaB],[],[x,y,psi0,L1]); 
    phi0=v(1);thetaB=v(2);
    
%nondimensional quantities
    FL2EIx=L1k^2*cos(thetaB); %nondimensional force Fx
    FL2EIy=L1k^2*sin(thetaB); %nondimensional force Fy
    MLEI=L3k*L1k; %nondimensional moment Frz
    
%output
    fx=-P*cos(thetaB);
    fy=-P*sin(thetaB);
    m=L3k*fy/k;
    x=L1minusXoverL1*L1-L1;
    y=-YoverL1*L1;
    rz=-psi0;
    F=P;
    p_=p;
    k_=k;
    phib_=phib;

%clear the global vars declared above
    clear L1k L3k L1minusXoverL1 YoverL1 k P E I FL2EIx FL2EIy MLEI p phib

%==========================

function closeness=nonlineardeflection(v,w)
    phi0=v(1);thetaB=v(2);
    x=w(1);y=w(2);psi0=w(3);L1=w(4);
global L1k L3k L1minusXoverL1 YoverL1 closeness k P E I FL2EIx FL2EIy MLEI p phib
%parameter definition
    theta0=thetaB+psi0; 
    p=sin(theta0/2);
    L3k=2*p*cos(phi0); %associated with moment 
    phib=asin(sin(thetaB/2)/p);
%elliptic integration
    Ephib=quadl(@ellipticintegrand,0,phib,[],[],p,2); %elliptic integral of the second kind 
    Fphib=quadl(@ellipticintegrand,0,phib,[],[],p,1); %elliptic integral of the fist kind 
    Ephi0=quadl(@ellipticintegrand,0,phi0,[],[],p,2); %elliptic integral of the second kind 
    Fphi0=quadl(@ellipticintegrand,0,phi0,[],[],p,1); %elliptic integral of the first kind 
    k=(Fphi0-Fphib)/L1; %associated with force
    L1k=L1*k;
    P=k*k*E*I;
    FL2EIx=L1k^2*cos(thetaB); %nondimensional force Fx
    FL2EIy=L1k^2*sin(thetaB); %nondimensional force Fy
    MLEI=L3k*L1k; %nondimensional moment Frz
    L1k=Fphi0-Fphib; %associated with force
    L1minusXoverL1=2*p*(cos(phib)-cos(phi0))/L1k; %beam shortening X
    YoverL1=(2*Ephib-Fphib)/L1k - (2*Ephi0-Fphi0)/L1k; %beam transverse deflection X
%error
    errormagnifier=1e3;
    closeness=( abs((L1minusXoverL1*L1-L1)+x) + abs(YoverL1*L1-y) )*errormagnifier;

%==========================

function F=ellipticintegrand(x,p,kind)
if kind==1 
   F=1./sqrt(1-(p.*sin(x)).^2); %integrand of the first kind
else
   F=sqrt(1-(p.*sin(x)).^2); %integrand of the second kind
end

%================================================================================================================================
%================================================================================================================================

% x_orig=-1e-6;y_orig=-5e-6;rz_orig=-pi/8 ;[phi0,thetaB,psi0,x,y,rz,fx,fy,m,closeness]=nonlinear_restoring_force2(x_orig,y_orig,rz_orig)
% x_orig=-0.5e-6;y_orig=-5e-6;rz_orig=-pi/8 ;[phi0,thetaB,psi0,x,y,rz,fx,fy,m,F,closeness]=nonlinear_restoring_force2(x_orig,y_orig,rz_orig),thetaB_pi2=thetaB/(pi/2),rz_pi2=rz/(pi/2),phi0_pi2=phi0/(pi/2)

% [x,y,rz,fx,fy,m,lm,p,k,phib]=NL_angles2Fq2(phi0,thetaB,psi0);dx=x,y,rz

% x_orig=-0.5e-6;y_orig=5e-6;rz_orig=-pi/8 ;[phi0,thetaB,psi0,x,y,rz,fx,fy,m,F,closeness]=nonlinear_restoring_force2(x_orig,y_orig,rz_orig),thetaB_pi2=thetaB/(pi/2),rz_pi2=rz/(pi/2),phi0_pi2=phi0/(pi/2),[x,y,rz,fx,fy,m,lm,p,k,phib]=NL_angles2Fq2(phi0,thetaB,rz);X=x,Y=y,RZ=rz
% p.fx=fx;p.fy=fy;p.m=m; net=cho_load('can06.m',p); q=cho_dc(net); figure(1); cho_display(net,q); qx=q(lookup_coord(net,'B','x')),qy=q(lookup_coord(net,'B','y')),qrz=q(lookup_coord(net,'B','rz'))
% [X,Y,RZ,S,DX,i]=NL_xypsi(p,k,phib,rz);figure(2);plot(X,Y,'g');

