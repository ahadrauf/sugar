%[x,y,psi0]=hjelm_vs_cho_Fy_2

%This plots displacement as a function of load, nondimensionally.

%Commands:
%X=0;Y=0;PSI=0;i=0;for p=0.1:0.1:25, i=i+1; [X(i),Y(i),PSI(i)]=hjelm_vs_cho_Fy_1(p); P(i)=p; end;
%figure(1);clf;plot(Y,P,'b',X,P,'r',PSI,P,'g');grid on;


function [X,Y,PSI]=hjelm_vs_cho_Fy_2(F)
tic;

nu = 0.3;
E = 165e9;
G = E / (2 * (1 + nu));
w = 2e-6;
h = 2e-6;
A = w*h;
I = h * w^3 / 12;
GA = G*A;
EI = E*I;
EA = E*A;
L = 20e-6;
%F = 1e-6;
%M = 1e-8*0;

%M = M;
Fy = F;
Fx = 0;
L1 = L;
EI = EI;

%Fy = P/L/L*EI;

%find p
ep = 0.001;
thetaB = pi/2;
psi01 = -pi/2+ep;
psi02 = pi/2-ep;
p01 = sin((thetaB + psi01)/2); %initial guess above
p02 = sin((thetaB + psi02)/2); %initial guess below

k = sqrt(Fy/E/I);
L1k = L1 * k;
p = fminbnd('elliptic_p_of_K_minus_Fphi',p01,p02,optimset('TolX',1e-15),L1k);
%p = fminbnd('elliptic_p_of_K_minus_Fphi',p01,p02,[],L1k);
phiB = asin(1/p/sqrt(2));
X = L1 - (2*p/k)*cos(phiB);
EphiB=quadl('elliptic',0,phiB,[],[],p,2);
FphiB=quadl('elliptic',0,phiB,[],[],p,1);
[Kp,Ep]=ellipke(p*p);
Y = ((2*EphiB - FphiB) - (2*Ep - Kp))/k;
PSI = 2*asin(p) - thetaB;

time=toc;

%normalize
X=X/L1;
Y=Y/L1;
PSI=PSI/(pi/2);
