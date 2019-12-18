function t=tofz(zoft,z,S,alpha,solution_guess)
%zoft='S*(t^(1-alpha)/(1-alpha)+t^(-alpha)/alpha)-z';
%S=3;alpha=0.7;z=4i;

t=solution_guess;
size_of_interval_guess=t*10;

dt=size_of_interval_guess;
dt0=dt;
dt1=dt0;
f0=eval(zoft);
I=0;

while abs(dt1)>1e-6
   t=t+dt0;
   f1=eval(zoft);
   dt1=-dt0*f1/(f1-f0);
   dt0=dt1;
   f0=f1;
   I=I+1;
end

iterations=I
