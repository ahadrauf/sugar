function sc14
resolution=50;
i=0;

g=2e-6;
pmin=0.19e-6;
pmax=1.4e-6;
dp=(pmax-pmin)/resolution;

for p=pmin:dp:pmax
    
    i=i+1;

    Zx1=p;
    Zy1=g+eps;
    t=-0.5;
    [t1,z1]=sc112(Zx1,Zy1,t,g,p);
    
    Zx2=p;
    Zy2=g+1e-6;
    [t2,z2]=sc112(Zx2,Zy2,t,g,p);
    e0=8.854e-12;
    h=40e-6;

   
    A=h*abs(Zy1-Zy2);
    Csc(i)=abs(e0*h/pi*log(t2/t1));
    Cpp(i)=A*e0/p;
    
    V=15;
    Fsc(i)=Csc(i)^2  * V^2 / (2*e0*A);
    Fpp(i)=Cpp(i)^2  * V^2 / (2*e0*A);
    Fpp(i)=1/2*A*e0* V^2 / p^2;

    x(i)=p;  
    
    Qsc(i)=Csc(i) * V;
    Qpp(i)=Cpp(i) * V;
end

figure(1);plot(x,2*Fpp,'g',x,2*Fsc,'b');title('Fringing field Force vs parallel comb, V=15, h=40um,w=2um');ylabel('Fsc (blue),  Fpp (green),  [N]');xlabel('g2 spacing,  [m]');grid on;
%figure(2);plot(x,-2*Fpp+2*Fsc,'r');title('fringing force minus parallel comb, V=15, h=40um,w=2um');ylabel('Fsc-Fpp (red)  [N]');xlabel('g2 spacing,  [m]');grid on;
figure(2);plot(x,Fsc./Fpp,'r');title('Fsc / Fpp, V=15, h=40um,w=2um');ylabel('Fsc/Fpp (red)  [N]');xlabel('g2 spacing,  [m]');grid on;

%figure(1);hold on;plot(rX,2*Fpp,'g',rX,2*Fsc,'b');

%figure(2);plot(X,Fpp,'g',X,Fsc,'b');title('Fringing field Force vs parallel plate');ylabel('Fsc (blue),  Fpp (green),  [N]');xlabel('length along plate,  [m]');grid on;
%figure(3);plot(X,Qpp,'g',X,Qsc,'b');title('Fringing field charge vs parallel plate');ylabel('Qsc (blue),  Qpp (green),  [C]');xlabel('length along plate,  [m]');grid on;

%figure(4);plot(0,0.011,'.w',1,t1,'o',2,t2,'o',3,t3,'o',4,t4,'o',5,-0.001,'.w');grid on;title('CPU times for total C, total F, total Q of fringing fields.');ylabel('time,  [seconds]');xlabel('Csc-blue, Fsc=green, Qsc=red, Cpp=cyan');
%figure(5);plot(0,0.011,'.w',1,t1,'o',2,t2,'o',3,t3,'o',4,t4,'o',6,-0.001,'.w',5,58,'o');grid on;title('CPU times for total C, total F, total Q of fringing fields.');ylabel('time,  [seconds]');xlabel('Csc-blue, Fsc=green, Qsc=red, Cpp=cyan, FEM=magenta');



%function sc11(gap,Zx,Zy,t)
%for gap closing actuator
%Zx and Zy are the physical coordinates
%The edge of the gap is on the right.
%t is the initial guess.

function [t,z]=sc112(Zx,Zy,t,g,p)
%By jvclark - May2002

    z=Zx + sqrt(-1)*Zy;
    [t]=fminsearch(@gap_type2,t,[],z,g,p);
    u=sqrt((t-(g/p)^2)/(t+1));
    z=2*g/pi*(p/g*atan(p*u/g) + 1/2*log((1+u)/(1-u)));

function [error_squared] = gap_type2(t,z,g,p)
    error_magnifier=1e30;
    u=sqrt((t-(g/p)^2)/(t+1));
    Z=2*g/pi*(p/g*atan(p*u/g) + 1/2*log((1+u)/(1-u)));
    error_squared = error_magnifier*(imag(z - Z))^2;

%[t]=fminsearch(   inline('(   (-1e-6 + sqrt(-1)*2e-6)     -(2e-6/pi)*(1+   (t)   +log(   (t)   )))^2'),  -0.5),  Z=(2e-6 / pi)*(1+t+log(t))
