%sc5

h=10e-6;
w=1e-6;
Csc=0;
Cpp=0;
G=0;
i=0;

g1=1e-6;
g2=6e-6;
L=5e-6;
I=sqrt(-1);
h=10e-6;
t1=-1;
E=8.854e-12;
for g=g1:(g2-g1)/15:g2   
   i=i+1;
   G(i)=g;
   Z=-L + I*g;
   T=sctdtfind(Z,w,g);
   t2=secant(Z,7*T,abs(T*5),'zoft_gap2',eps,g,w)
   Csc(i)=2*h*E/pi*log(t1/t2);
   Cpp(i)=2*E*L*h/g;
end

figure(1);plot(G,real(Csc),'b',G,real(Cpp),'r');ylabel('CFringe-blue  CApprox-reb  [F/m]');xlabel('gap [meters]');
figure(2); plot(G,real(((Csc)-(Cpp))/(Csc)),'b');

