%[x,y,psi0]=hjelm_vs_cho1
%hjelm_vs_cho1

function [x,y,psi0]=hjelm_vs_cho1
tic;
M = 0;
Fy = -5;
Fx = 5;
F0 = sqrt(Fx^2 + Fy^2);
L3 = M/F0;
L1 = 10;
EI = 1000;
thetaB = acos(Fx/F0);
thetaB = asin(Fy/F0);
k = sqrt(F0/EI);
L3k = L3*k;
L1k = L1*k;
thetaB2 = thetaB/2;

%find p
ep = 0.001;
psi01 = -pi/2+ep;
psi02 = pi/2-ep;
p01 = sin((thetaB + psi01)/2); %initial guess above
p02 = sin((thetaB + psi02)/2); %initial guess below

p = fminbnd('elliptic_p_of_Fphi0_minus_FphiB',p01,p02,[],thetaB2,L3k,L1k);
phiB = asin(sin(thetaB/2)/p);
phi0 = acos(L3k/2/p);
x = L1 - 2*p*(cos(phiB) - cos(phi0))/L1k*L1;
EphiB=quadl('elliptic',0,phiB,[],[],p,2);
FphiB=quadl('elliptic',0,phiB,[],[],p,1);
Ephi0=quadl('elliptic',0,phi0,[],[],p,2);
Fphi0=quadl('elliptic',0,phi0,[],[],p,1);
y = L1*((2*EphiB - FphiB) - (2*Ephi0 - Fphi0))/L1k;
psi0 = 2*asin(p) - thetaB;
time=toc
