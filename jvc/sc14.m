function sc14
resolution=100;
xmax=1e-6;
dx=xmax/resolution;
i=0;


gmin=2e-6;
gmax=0.2e-6;
dg=(gmax-gmin)/resolution;

for g=gmin:dg:gmax
    
    i=i+1;
    %g=2e-6;

    Zx1=-1e-6;
    Zy=g;
    t=-0.5+eps;
    [t1,z1]=sc11(g,Zx1,Zy,t);

    Zx2=-eps;
    [t2,z2]=sc11(g,Zx2,Zy,t);
    
    e0=8.854e-12;
    h=20e-6;

    A=h*abs(Zx1-Zx2);
    Csc(i)=abs(e0*h/pi*log(t2/t1));
    Cpp(i)=A*e0/g;
    
    V=15;
    Fsc(i)=Csc(i)^2  * V^2 / (2*e0*A);
    Fpp(i)=Cpp(i)^2  * V^2 / (2*e0*A);
    Fpp(i)=1/2*A*e0* V^2 / g^2;

    x(i)=g;  
    
    Qsc(i)=Csc(i) * V;
    Qpp(i)=Cpp(i) * V;
end

figure(1);plot(x,2*Fpp,'g',x,2*Fsc,'b');title('Fringing field Force vs parallel plate, V=15, L=20um');ylabel('Fsc (blue),  Fpp (green),  [N]');xlabel('gap spacing,  [m]');grid on;
%figure(2);plot(x,-2*Fpp+2*Fsc,'r');title('fringing force minus parallel plate, V=15, L=20um, w=h');ylabel('Fsc-Fpp (red)  [N]');xlabel('gap spacing,  [m]');grid on;
figure(2);plot(x,Fsc./Fpp,'r');title('Fsc/Fpp, V=15, L=20um, w=h');ylabel('Fsc/Fpp (red)  [N]');xlabel('gap spacing,  [m]');grid on;
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

function [t,z]=sc11(g,Zx,Zy,t)
%By jvclark - May2002

%eg: sc11(2e-6,-1e-6,2e-6,-0.5)

S=g/pi;

t=-0.5;
z=Zx + sqrt(-1)*Zy;
[t]=fminsearch(@gap_type1,t,[],z,S);
z=S*(1+t+log(t));

function [error_squared] = gap_type1(t,z,S)
    error_magnifier=1e20;
    error_squared = error_magnifier*(z - S*(1 + t + log(t)))^2;

%[t]=fminsearch(   inline('(   (-1e-6 + sqrt(-1)*2e-6)     -(2e-6/pi)*(1+   (t)   +log(   (t)   )))^2'),  -0.5),  Z=(2e-6 / pi)*(1+t+log(t))
