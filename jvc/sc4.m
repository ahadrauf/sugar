
d1=0.01e-6;
d2=10e-6;
step=(d2-d1)/100
i=0;
t2=-1;
for d=d1:step:d2
   i=i+1;
   z1(i)=eval('d*(1+t+log(t))/pi');
   Csc(i)=h*8.854e-12/pi*log(-1/t(i));
   Cpp(i)=h*8.854e-12*real(z)/d;


