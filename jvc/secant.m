function r=secant(z,t_guess,t_range,fnc,small,gap,width)
%gap=0.9e-6;t2=secant(-10e-6+I*2e-6,-1e-8,1e-10,'zoft_gap2',eps,gap);t1=-1;Csc=8.854e-12/pi*log(t1/t2),Cpp=8.854e-12*10e-6/gap
%t=secant(50e-6+I*2e-6,-1e-10,0.01,'zoft_gap2',eps)
%i=0;for g=1e-6:(-1e-6 + 6e-6)/100:6e-6;i=i+1;G(i)=g;gap=4e-6;t2=secant(-10e-6+I*2e-6,-1e-10,1e-15,'zoft_gap2',eps,G(i),G(i));t1=-1;Csc(i)=8.854e-12/pi*log(t1/t2);Cpp(i)=8.854e-12*10e-6/G(i);end
% figure(1);plot(G,Cpp);

%width=1e-6;
%gap=2e-6;
r=t_guess;
dx0=t_range;
dx1=dx0;
f0=feval(fnc,width,gap,r);
f0=f0-z;
i=0;

while abs(dx1)>small
   r=r+dx0;
   f1=feval(fnc,width,gap,r);
   f1=f1-z;
   dx1=-dx0*f1/(f1-f0);
   if abs(dx1)>1e4; error('method diverging');end
   dx0=dx1;
   f0=f1;
   i=i+1;
   if i>2000; error('TOO MANY INTERATIONS');end
end
iterations=i;

