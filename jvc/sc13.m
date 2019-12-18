function sc13
resolution=100;
xmax=1e-6;
dx=xmax/resolution;
i=0;

for x=0:dx:xmax-dx

    i=i+1;
    g=2e-6;

    Zx1=-(xmax-x + eps);
    Zy=g;
    t=-10+eps;
    tic;
    [t1,z1]=sc11(g,Zx1,Zy,t);

    Zx2=-(xmax-x-dx+eps);
    [t2,z2]=sc11(g,Zx2,Zy,t);
    
    e0=8.854e-12;
    h=20e-6;

    A=h*abs(Zx1-Zx2);
    Csc(i)=abs(e0*h/pi*log(t2/t1));
    t1=toc;
    tic;
    Cpp(i)=A*e0/g;
    t2=toc;
    X(i)=x;
    rX(resolution-i+1)=x;
    
    V=15;
    tic;
    Fsc(i)=Csc(i)^2  * V^2 / (2*e0*A);
    t3=toc;
    Fpp(i)=Cpp(i)^2  * V^2 / (2*e0*A);
    
    tic
    Qsc(i)=Csc(i) * V;
    t4=toc;
    Qpp(i)=Cpp(i) * V;
end

figure(1);plot(xmax+X,Cpp,'g',xmax+X,Csc,'b');title('Fringing Force vs parallel plate, V=15,L=20um,w=h=2um=g');ylabel('Fsc (blue),  Fpp (green),  [N]');xlabel('width along gap,  [m]');grid on;
figure(1);hold on;plot(rX,Cpp,'g',rX,Csc,'b');

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
