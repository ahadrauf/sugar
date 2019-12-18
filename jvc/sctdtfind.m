function [t]=sctdtfind(z,w,g)
width=w;
gap=g;

t=-0.5;
p=width;
g=gap;
k=1+p/g;
a=-1+2*k^2+ 2*k*sqrt(k^2-1);
R=sqrt((t+1)/(t+a));
forward=1;

while forward>0
   R=sqrt((t+1)/(t+a));   
   zt=g/pi*((a+1)/sqrt(a)*atanh(R)+(a-1)/sqrt(a)*R/(1-R^2)+log((R*sqrt(a)-1)/(R*sqrt(a)+1)));
   if abs(real(zt))>abs(real(z))
      forward=0;
   end
   t=t/2;
end

%Csc=0;Cpp=0;G=0;i=0;for g=1e-6:(-1e-6 + 5e-6)/10:5e-6, i=i+1;G(i)=g;Z=-5e-6 + I*g;T=sctdtfind(Z,G(i));t2=secant(Z,T,T*5,'zoft_gap2',eps,G(i),2e-6); t1=-1;Csc(i)=8.854e-12/pi*log(t1/t2);Cpp(i)=8.854e-12*10e-6/G(i);end;

%h=10e-6;Csc=0;Cpp=0;G=0;i=0;for g=1e-6:(-1e-6 + 5e-6)/10:5e-6, i=i+1;G(i)=g,;Z=-5e-6 + I*g;T=sctdtfind(Z,G(i));t2=secant(Z,T,T*5,'zoft_gap2',eps,G(i),2e-6); t1=-1;Csc(i)=h*8.854e-12/pi*log(t1/t2);Cpp(i)=8.854e-12*5e-6*h/G(i);end;
