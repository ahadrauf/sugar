function t=tofz(z,t_guess,t_range)
%(zoft,z,S,alpha,solution_guess)
%zoft='S*(t^(1-alpha)/(1-alpha)+t^(-alpha)/alpha)-z';
%S=3;alpha=0.7;z=4i;

d=2e-6;
t=t_guess;
%Z=eval('d*(1+t+log(t))/pi')
zoft='d*(1+t+log(t))/pi-z';
%zoft='d*(1+t+log(t))/pi';
%z=eval(zoft)
%z=0+i*d;
%d=2;

dt=t_range;
dt0=dt;
dt1=dt0;
f0=eval(zoft);
I=0;
tol=eps*2;
while abs(dt1)>tol
   t=t+dt0;
   f1=eval(zoft);
   dt1=-dt0*f1/(f1-f0);
   dt0=dt1;
   f0=f1;
   I=I+1;
end
iterations=I

TT=t
ZZ=eval(zoft)