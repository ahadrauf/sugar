%function [minimization] = minimize_fn_theta1_phi0(theta1,phi0,XL1,YL1,psi)
%For a nonlinear beam subject to simultaneous Fx, Fy, M
%For use with fmins to find theta1 and phi0 such that Y()=Y, X()=X, psi()=psi.
%theta1 is the angle of the external node force.
%phi0 is the projection angle of the moment arm.
%XL1 is the projected beam shortening per beam length.
%YL1 is the lateral deflection per beam length.
%psi is the angular displacement.
%See also elliptic,

%op=foptions;op(2)=1e-2;op(3)=1e-2;
%L1=100e-6; theta1=pi/4; phi0=pi/4; XL1=2e-6/L1; YL1=20e-5/L1;psi=pi/5;
%thetaphi=fmins('minimize_fn_theta1_phi0',[theta1,phi0],options,[],XL1,YL1,psi)
 
function [minimization] = minimize_fn_theta1_phi0(parameter,XL1,YL1,psi)
%By: JVClark - Dec2001.

%Uses: elliptic,FxFyM_fn_theta1_phi0
%Used by:

theta1=parameter(1);
phi0=parameter(2);

theta0=theta1+psi;
p=sin(theta0/2); %moduls of the elliptic integral.
phi1=asin(sin(theta1/2)/p); 

Ephi0=quad8('elliptic',0,phi0,[],[],p,2); %elliptic integral of the 2nd kind
Ephi1=quad8('elliptic',0,phi1,[],[],p,2);
Fphi0=quad8('elliptic',0,phi0,[],[],p,1); %elliptic integral of the 1st kind
Fphi1=quad8('elliptic',0,phi1,[],[],p,1);

L1k=Fphi0-Fphi1;

XL1_of_theta1_phi0=2*p*(cos(phi1)-cos(phi0))/L1k;
YL1_of_theta1_phi0=(2*Ephi1-Fphi1)/L1k - (2*Ephi0-Fphi0)/L1k;
psi_of_theta1_phi0=theta0-theta1;

minimization=abs(abs(XL1_of_theta1_phi0-XL1)+(YL1_of_theta1_phi0-YL1)+(psi_of_theta1_phi0-psi));

