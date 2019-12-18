function [t,z]=sc_parallel_newton(z,t)
small=eps;
small=1e-50;
I=sqrt(-1);
dt=1;
i=0;

d=1e-6;
zoft='d*(1+t+log(t))/pi-z';
dzdt='d*(1+1/t)/pi';

while abs(dt) > small
   i=i+1;   
   if i>20000, error(' ::: EXCESSIVE ITERATIONS ::: '); end
   df=eval(dzdt);
   if abs(df) < small, error(' ::: DIVERGING ::: '); end
   dt=-eval(zoft)/df;
   t=t+dt;
end
iterations=i

t
ZZ=eval('d*(1+t+log(t))/pi')