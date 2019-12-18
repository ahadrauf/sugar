%function Z=nonlineardeflection(x,y,rz)
%Given x, y, and rz displacements, 
%determine angles thetaB and phi0 (see derivation for term identification).

function [FL2EIx,FL2EIy,MLEI]=nonlinear_restoring_force(x,y,rz)
%jvclark - April 2002 

%global vars ONLY for the functions below
global L1k L3k L1minusXoverL1 YoverL1 Z

%starting guess
phi0 = pi/2 - pi/10;
thetaB = pi/2 - pi/10;
psi0 = rz;

phi0 = pi/2 - pi/10;
thetaB = pi/2 - pi/10;
psi0 = rz;

%true angles of forces
v=fminsearch(@nonlineardeflection,[phi0,thetaB,psi0],[],[x,y]);
phi0=v(1);thetaB=v(2);psi0=v(3);

%forces
FL2EIx=L1k^2*cos(thetaB); %nondimensional force Fx
FL2EIy=L1k^2*sin(thetaB); %nondimensional force Fy
MLEI=L3k*L1k; %nondimensional moment Frz


X=-(L1minusXoverL1*L1k-L1k)
Y=YoverL1*L1k 
Z

%clear the global vars declared above
clear L1k L3k L1minusXoverL1 YoverL1 Z


%==========================

function F=elliptic(x,p,kind)
if kind==1 
   F=1./sqrt(1-(p.*sin(x)).^2); %integrand of the first kind
else
   F=sqrt(1-(p.*sin(x)).^2); %integrand of the second kind
end

%============================

function Z=nonlineardeflection(v,w)
phi0=v(1);thetaB=v(2);psi0=v(3);
x=w(1);y=w(2);
global L1k L3k L1minusXoverL1 YoverL1 Z
%parameter definition
theta0=thetaB+psi0; 
p=sin(theta0/2);
L3k=2*p*cos(phi0); %associated with moment 
phib=asin(sin(thetaB/2)/p);
%elliptic integration
Ephib=quadl('elliptic',0,phib,[],[],p,2); %elliptic integral of the second kind 
Fphib=quadl('elliptic',0,phib,[],[],p,1); %elliptic integral of the fist kind 
Ephi0=quadl('elliptic',0,phi0,[],[],p,2); %elliptic integral of the second kind 
Fphi0=quadl('elliptic',0,phi0,[],[],p,1); %elliptic integral of the first kind 
%parameter definition
L1k=Fphi0-Fphib; %associated with force
L1minusXoverL1=2*p*(cos(phib)-cos(phi0))/L1k; %beam shortening X
YoverL1=(2*Ephib-Fphib)/L1k - (2*Ephi0-Fphi0)/L1k; %beam transverse deflection X
Z=abs(-(L1minusXoverL1*L1k-L1k)-x)+abs(YoverL1*L1k-y)*1000000;

%================================================================================================================================
%================================================================================================================================
%================================================================================================================================

