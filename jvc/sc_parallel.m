function [F,C,t1,t2]=sc_parallel(z1,z2)

t=tofz(zoft,z,S,alpha,solution_guess);

t=solution_guess;
size_of_interval_guess=t*10;

dt=size_of_interval_guess;
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
