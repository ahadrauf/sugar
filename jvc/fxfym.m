% [Y,X,Fx,Fy,M]=fxfym(pi/2,pi/2,psi0)  %thetab,phi0,psi0
function [YoverL1,L1minusXoverL1,FL2EIx,FL2EIy,MLEI]=fxfym(thetab,phi0,psi0) 
theta0=thetab+psi0; 
p=sin(theta0/2); %parameter 
L3k=2*p*cos(phi0); %associated with moment 
phib=asin(sin(thetab/2)/p);

%elliptic integration
Ephib=quad8('elliptic',0,phib,[],[],p,2); %elliptic integral of the second kind 
Fphib=quad8('elliptic',0,phib,[],[],p,1); %elliptic integral of the fist kind 
Ephi0=quad8('elliptic',0,phi0,[],[],p,2); %elliptic integral of the second kind 
Fphi0=quad8('elliptic',0,phi0,[],[],p,1); %elliptic integral of the first kind 

L1k=Fphi0-Fphib; %associated with force
L1minusXoverL1=2*p*(cos(phib)-cos(phi0))/L1k; %beam shortening 
YoverL1=(2*Ephib-Fphib)/L1k - (2*Ephi0-Fphi0)/L1k; %beam transverse deflection 
FL2EIx=L1k^2*cos(thetab); %nondimensional force
FL2EIy=L1k^2*sin(thetab); %nondimensional force
MLEI=L3k*L1k; %nondimensional moment

