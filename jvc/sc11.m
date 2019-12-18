%function sc11(gap,Zx,Zy,t)
%for gap closing actuator
%Zx and Zy are the physical coordinates
%The edge of the gap is on the right.
%t is the initial guess.

function [t,z]=sc11(g,Zx,Zy,t)
%By jvclark - May2002

%eg: sc11(2e-6,-1e-6,2e-6,-0.5)

S=g/pi;

t=-0.5;
z=Zx + sqrt(-1)*Zy;
[t]=fminsearch(@gap_type1,t,[],z,S);
z=S*(1+t+log(t));

function [error_squared] = gap_type1(t,z,S)
    error_magnifier=1e10;
    error_squared = error_magnifier*(z - S*(1 + t + log(t)))^2;

%[t]=fminsearch(   inline('(   (-1e-6 + sqrt(-1)*2e-6)     -(2e-6/pi)*(1+   (t)   +log(   (t)   )))^2'),  -0.5),  Z=(2e-6 / pi)*(1+t+log(t))
