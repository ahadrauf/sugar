clear all
p=2e-6; 
g=2e-6; 
V=50; 
h=50e-6;   
epsilon=8.854e-12;
%z1=p*7+j*g; 
%z2=p*8+j*g;  

i=0;
P=p+eps+0e-6;
dp=1.5e-6;
for g=6e-6 : -0.1e-6 : 1.1e-6
   i=i+1;
   G(i)=g;
   
   z1=p+P    +j*(g);
   z2=p+P+dp +j*(g);
   
   [C(i),Q(i),F(i),Ca(i),Fa(i)]=sc_force(V,h,p,g,[z1,z2],dp);
   
end
figure(2);
clf;
plot(G,F*2,'b');xlabel('gap2 [meter]');ylabel('Force [Newtons]');
hold on;
plot(G,Fa,'r');

figure(1);
clf;
plot(G,C*2,'b');xlabel('gap2 [meter]');ylabel('Capacitance [Farads/meter]');
hold on;
plot(G,Ca,'r');
