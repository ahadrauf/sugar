clear all;
%plot1
%g1=0.6e-6;
%g2=1e-6;
%plot2
%g1=1e-6;
%g2=10e-6;
resolution=40;
W=1e-6;
H=10e-6;
L1=0e-6;
L2=20e-6;
L=L2-L1;
tol=1e-20;
E=8.845e-12;
t1=-1;
i=0;
V=50;
A=L*H;

for L1 = 0:5e-6 / 80:5e-6
   L2 = L1 + 5e-6 / 80;
   g=1e-6;  
   L=L2-L1;
   
%for g=g1:(g2-g1)/resolution:g2
   i=i+1;
   z1=-L1 + sqrt(-1)*g;
   z2=-L2 + sqrt(-1)*g;
   z=-L + sqrt(-1)*g;
   t2=sc7(W,g,z1,tol);
   t1=sc7(W,g,z2,tol);
   G(i)=g;
   Cpp(i)=2*E*H*L/(g);
   Csc(i)=2*H*E/pi*log(t2/t1);
   C=(Csc(i)-Cpp(i))*L/H;
   Csc(i)=Csc(i)+C;
   Qsc=Csc(i)*V;
   Qpp=Cpp(i)*V;
   Fsc(i)=2 * Qsc^2/(4*pi*E*g^2);
   Fpp2(i)=2 * Qpp^2/(4*pi*E*g^2);
   Fpp(i)=2 * 0.5*E*A*V^2/g^2;
   LL(i)=L1;
end

%plots 1 and 2
%figure(1);plot(G,Cpp,'-b',G,Csc,'r');xlabel('gap [m]');ylabel('C(Fringing)-red  C(ParallelApprox)-blue  [F/m]');title('C vs gap for parallel fingers');
%figure(2);plot(G,(Csc-Cpp)./Csc,'k');xlabel('gap [m]');ylabel('[C(Fringing)-C(ParallelApprox)]/C(Fringing)  [unitless]');title('(Cf-Cp)/Cf vs gap for parallel fingers');
%figure(3);plot(G,Fpp,'b',G,Fsc,'r',G,Fpp2,'g');xlabel('gap [m]');ylabel('F(Fringing)-red  F(ParallelApprox)-blue  [F/m]');title('F vs gap for parallel fingers');

%plot3
figure(2);plot(LL,Cpp/L,'b',LL,Csc/L,'r');xlabel('x along beam length [m]');ylabel('CFringe/m-red  CApprox/m-blue  [F/m/m]');title('C per unit length. Parallel electrodes');
