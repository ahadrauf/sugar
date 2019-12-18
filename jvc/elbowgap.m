function [f]=elbowgap(t,z)

t
z

g=2e-6;
p=4e-6;

n = t - (g/p)^2;
d = t + 1;
u = sqrt(n/d);

f = -z + 2*g/pi*( p/g*atan(p*u/g) + 1/2*log((1+u)/(1-u)));

f