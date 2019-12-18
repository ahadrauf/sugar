clear all;
E=160e9;
nu=0.226;
G=E/(2*(1+nu));
L=50e-6;
h=2e-6;
w=2e-6;
I=w^3*h/12;
L3EI3=L^3/(3*E*I);
F=30e-6; 

%rect
i=0;
i=i+1;
y1(i)=6.04665e-6; r1(i)=32e-6;
i=i+1;
y1(i)=6.04675e-6; r1(i)=16e-6;
i=i+1;
y1(i)=6.04587e-6; r1(i)=8e-6;
i=i+1;
y1(i)=6.03979e-6; r1(i)=4e-6;
i=i+1;
y1(i)=6.00969e-6; r1(i)=2e-6;

%with 1 fixed boundary
%i=i+1;
%y1(i)=6.19987e-6; r1(i)=1e-6; 
Dy50=30e-6*(50e-6)^3 / (3*E*I) * 1.0029
i=i+1;
r1(i)=0e-6;
y1(i)=Dy50;

%circ
i=0;
i=i+1;
y2(i)=5.99650e-6; r2(i)=1e-6;
i=i+1;
y2(i)=6.03397e-6; r2(i)=2e-6;
i=i+1;
y2(i)=6.06054e-6; r2(i)=4e-6;
i=i+1;
y2(i)=6.07040e-6; r2(i)=8e-6;
i=i+1;
y2(i)=6.10109e-6; r2(i)=16e-6;
%half circle
%y2(i)=6.24526e-6; r(i)=10e-6; 

%plot(r2,y2,'-*r');
grid on;

%--------------------------------------------------------------------

E=160e9;
nu=0.226;
G=E/(2*(1+nu));
L=50e-6;
h=2e-6;
w=2e-6;
I=w^3*h/12;
L3EI3=L^3/(3*E*I);
F=30e-6; 
T=L*F; 
d=2e-6;


T=T*2*2;

i=0;
i=i+1;
r(i)=32e-6;
D=r(i)*2;
phi=T*((1/d^2 - 1/D^2)/(pi*h*G)); 
Y(i)=F*L3EI3+L*sin(phi);

i=i+1;
r(i)=16e-6;
D=r(i)*2;
phi=T*((1/d^2 - 1/D^2)/(pi*h*G)); 
Y(i)=F*L3EI3+L*sin(phi);

i=i+1;
r(i)=8e-6;
D=r(i)*2;
phi=T*((1/d^2 - 1/D^2)/(pi*h*G)); 
Y(i)=F*L3EI3+L*sin(phi);

i=i+1;
r(i)=4e-6;
D=r(i)*2;
phi=T*((1/d^2 - 1/D^2)/(pi*h*G)); 
Y(i)=F*L3EI3+L*sin(phi);

i=i+1;
r(i)=2e-6;
D=r(i)*2;
phi=T*((1/d^2 - 1/D^2)/(pi*h*G)); 
Y(i)=F*L3EI3+L*sin(phi);

i=i+1;
r(i)=1e-6;
D=r(i)*2;
phi=T*((1/d^2 - 1/D^2)/(pi*h*G))

r(i)=0;
Y(i)=F*L3EI3+L*sin(phi)

%i=i+1;
%r(i)=0e-6;
%D=r(i)*2;
%phi=T*((1/d^2 - 1/D^2)/(pi*h*G)); 
%Y(i)=F*(51e-6)^3 / (3*E*I) +L*sin(phi);

figure(1); clf;
plot(r1,y1,'-*b');
hold on;
plot(r,Y,'-*g');
grid on;
title('tip deflection vs base radius');
xlabel('base radius [m]    E=160e9,nu=0.226,w=h=2um,L=50um');
ylabel('tip deflection = Y  [m]    tip F=30uN');

figure(2);clf;
plot(r,(Y-y1)./Y,'-r*');
grid on;
title('relative error vs base radius');
xlabel('base radius [m]    E=160e9,nu=0.226,w=h=2um,L=50um');
ylabel('relative error = (Y-y)/Y');

figure(3);clf;
plot(r,(Y-y1),'-r*');
grid on;
title('error vs base radius');
xlabel('base radius [m]    E=160e9,nu=0.226,w=h=2um,L=50um');
ylabel('error = Y-y  [m]');

figure(4);clf;
plot(r(2:6),(Y(2:6)-Y(1:5))./Y(2:6),'-b*');
grid on; hold on;
plot(r(2:6),(y1(2:6)-y1(1:5))./y1(2:6),'-g*');
grid on;
title('relative convergence vs base radius');
xlabel('base radius [m]    E=160e9,nu=0.226,w=h=2um,L=50um');
ylabel('relative convergence = (y(2:6)-y(1:5))./y(2:6)  [m]');
