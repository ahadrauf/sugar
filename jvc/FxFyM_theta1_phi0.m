%function [FL2EIX,FL2EIY,MLEI] = FxFyM_fn_theta1_phi0(theta1,phi0,psi)
%For a nonlinear beam subject to simultaneous Fx, Fy, M
%Assuming you found theta1 and phi0 such that Y()=Y, X()=X, psi()=psi.
%theta1 is the angle of the external node force.
%phi0 is the projection angle of the moment arm.
%See also: minimize_fn_theta1_phi0

% L1=100e-6;theta1=pi/4;phi0=pi/4;XL1=2e-6/L1;YL1=20e-5/L1;psi=pi/5;[minimization] = minimize_fn_theta1_phi0(theta1,phi0,XL1,YL1,psi)
 
function [FL2EIX,FL2EIY,MLEI]=FxFyM_theta1_phi0(theta1,phi0,psi)
%By: JVClark - Dec2001.

%Uses: 
%Used by:

theta0=theta1+psi;
p=sin(theta0/2); %moduls of the elliptic integral.
phi1=asin(sin(theta1/2)/p); 

Fphi0=quad8('elliptic',0,phi0,[],[],p,1); %elliptic integral of the 1st kind
Fphi1=quad8('elliptic',0,phi1,[],[],p,1);

L3k=2*p*cos(phi0);
L1k=Fphi0-Fphi1;

FL2EIX=L1k^2*cos(theta1);
FL2EIY=L1k^2*sin(theta1);
MLEI=L3k*L1k;

