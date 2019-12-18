%p=2e-6; g=2e-6; V=50; h=50e-6;   z1=p*7+j*g; z2=p*8+j*g;  sc_force(V,h,p,g,[z1,z2])

function [C,Q,F,Capprox,Fapprox]=sc_force(V,h,p,g,z,dp)  
%By: Jason Vaughn Clark - Nov2001 
%uses sc_comb1

%t1=lower integral bound
%t2=upper intergral bound
%V=voltage
%h=layer thickness, out of the Schwarz-Christofel plane

epsilon0=8.854e-12; %[F/m] permitivity of free space

%Find the bounds on T
for i=1:2
   if real(z(i))>p %horizontal section of the moveable finger
      t_range=[-1-eps,-1e15]; %t range
      [T(i),Zero(i)]=fzero('sc_comb1',t_range,optimset('Display','off'),p,g,z(i)); %Find corresponding t for z.
   else %vertical section
      t_range=[-1+eps,1e15]; %t range
      [T(i),Zero(i)]=fzero('sc_comb1',t_range,optimset('Display','off'),p,g,z(i)); %Find corresponding t for z.
   end
end

%Calculate C capacitance, F force, and Q charge.
C=h*epsilon0/pi*log(T(2)/T(1)); %capacitance
Q=V*C; %charge
A=h*dp; %area
F=(Q/A)^2 * A/epsilon0/2; %force
C=C;
F=F;

A=h*dp;
Capprox=epsilon0*A/g;
Fapprox=0.5*epsilon0*A*V^2/g^2;

%
% |      | t=0
% |      |
% |      |
% |      |
% |  p   |
% |      |
% |      |
% |      | t=-1                    t=-inf
% |      +-------------------------------
% |
% |            jg
% |z=0
% +--------------------------------------
%                                  t=+inf

