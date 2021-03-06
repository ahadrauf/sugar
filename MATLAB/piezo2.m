E=4e9; cm=0.01;
h=2e-6; 
w=2e-6; 
L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; k2=k1/2; NH=cm/w, NW=cm/(2*L), ND=cm/h, K=k2*NH*NW*ND; 
Y=2.63e-5; 
m=1800*cm^3; h=K*Y^2/(2*9.8*m)
% L   w   y
%250 100 0.0192
%250  50 0.0114
%500 100 0.0164
%250 200 0.0305

%500 100 0.0047
%500 200 0.0030
%500  50 0.0040
%50    5 0.0040
%200   5 0.0040

%100   5 0.4848
%100  10 0.5895
%100  20 0.3153
%100  10 1.8980

%inc w, decL 

E=4e9; cm=0.01;
h=2e-6; 
w=10e-6; 
L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+5e-6), NW=cm/(4*L), ND=cm/h, k2=NW*k1; k3=k2/NH; k4=k3*ND; 
Y=1.846e-5; 
m=1800*cm^3; h=k4*(Y*NH)^2/(2*9.8*m)
%h=1.6
E=4e9; cm=0.01;
h=2e-6; 
w=10e-6; 
L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+5e-6), NW=cm/(4*L), ND=cm/h, k2=k1*2; k3=k2/NH; k4=k3*NW; K=k4*ND; 
Y=1.846e-5; 
m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m)
%h=3.2

E=4e9; cm=0.01;
h=2e-6; 
w=20e-6; 
L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+5e-6), NW=cm/(4*L), ND=cm/h, k2=k1*2; k3=k2/NH; k4=k3*NW; K=k4*ND; 
Y=9.783e-6; 
m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m)
%h=4.3

%E=4e9; h=6e-6; w=20e-6; A=w*h; L = 204e-6; x = (2e-6)/(E*A/L)

E=4e9; cm=0.01; h=6e-6; w=20e-6; L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+5e-6), NW=cm/(4*L), ND=cm/h, k2=k1*2; k3=k2/NH; k4=k3*NW; K=k4*ND; 
Y=9.07e-7; 
m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m)
%h=0.0373
E=4e9; cm=0.01; h=6e-6; w=20e-6; L=200e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+5e-6), NW=cm/(4*L), ND=cm/h, k2=k1*2; k3=k2/NH; k4=k3*NW; K=k4*ND; 
Y=3.576e-6;
m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m)
%h=0.0362
E=4e9; cm=0.01; h=6e-6; w=40e-6; L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+5e-6), NW=cm/(4*L), ND=cm/h, k2=k1*2; k3=k2/NH; k4=k3*NW; K=k4*ND; 
Y=5.858e-7;
m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m)
%h=0.0692

%E=4e9; h=2e-6; w=20e-6; A=w*h; L = 204e-6; x = (2e-6)/(E*A/L), %2.5500e-009
E=4e9; cm=0.01; h=2e-6; w=20e-6; L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+5e-6), NW=cm/(4*L), ND=cm/h, k2=k1*2; k3=k2/NH; k4=k3*NW; K=k4*ND; 
Y=1.165e-6;
m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m)
%h=0.0616

E=4e9; cm=0.01; h=2e-6; w=4e-6; L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+2e-6), NW=cm/(4*L), ND=cm/h, k2=k1*2; k3=k2/NH; k4=k3*NW; K=k4*ND; 
Y=2.453e-6;
m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m)
%h=0.0091

%h=2e-6; w=4e-6; A=w*h; L = 208e-6; x=1.34e-9;  E = (2e-6)/(x*A/L)
E=3.8806e10; cm=0.01; h=2e-6; w=4e-6; L=100e-6; 
I=w^3*h/12; k1=3*E*I/L^3; NH=cm/(w+2e-6), NW=cm/(4*L), ND=cm/h, k2=k1*2; k3=k2/NH; k4=k3*NW; K=k4*ND; 
Y=2.58e-6;
m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m)
%h=0.0976

L=200e-6; w=4e-6; h=2e-6;
K2=0.239; NW=cm/L; NH=cm/(w+2e-6); ND=cm/h; KW=NW*K2; KH=KW/NH; K=ND*KH;
Y=3.549e-6; m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m), J=0.5*K*(Y*NH)^2
C=8.854e-12 * 2e-6 * 200e-6 / 2e-6; V=5.4;
JJ=0.5*C*V^2
%0.0356

L=200e-6; w=40e-6; h=2e-6;
K2=219.9736; NW=cm/L; NH=cm/(w+2e-6); ND=cm/h; KW=NW*K2; KH=KW/NH; K=ND*KH;
Y=0.61e-6; m=1800*cm^3; h=K*(Y*NH)^2/(2*9.8*m), J=0.5*K*(Y*NH)^2
C=8.854e-12 * 20e-6 * 200e-6 / 2e-6; V=5.4;
JJ=0.5*C*V^2
%0.1381

'pvdf'
x=6.369e-5; f=1e-6; k = f/x; x=1.124e-5; J1 = 0.5 * k * x^2 / (200e-4 * 2e-4 * 4e-4);
Nm = k*x * x; L=200e-6; w=4e-6; h=2e-6; NW=cm/L; NH=cm/(w+2e-6); ND=cm/h; k=NW*k; k=k/NH; k=ND*k; m=1800*cm^3; 
H=k*(x*NH)^2/(2*9.8*m), J=0.5*k*(x*NH)^2

%pvdf
x=6.62e-8; f=1e-6; k = f/x;
x=1.124e-6; J1 = 0.5 * k * x^2 / (200e-4 * 2e-4 * 40e-4);
Nm = k*x * x;
L=200e-6; w=40e-6; h=2e-6;
NW=cm/L; NH=cm/(w+2e-6); ND=cm/h; k=NW*k; k=k/NH; k=ND*k;
m=1800*cm^3; H=k*(x*NH)^2/(2*9.8*m), 
J=0.5*k*(x*NH)^2
C=8.854e-12 * 40e-6 * 200e-6 / 2e-6; V=10; JJ=0.5*C*V^2;

'pzt'
x=5.11e-7; f=1e-6; k = f/x,
x=11.25e-6; J1 = 0.5 * k * x^2 / (200e-4 * 2e-4 * 4e-4);
Nm = k*x * x, L=200e-6; w=4e-6; h=2e-6; NW=cm/L; NH=cm/(w+2e-6); ND=cm/h; k=NW*k; k=k/NH; k=ND*k; m=7500*cm^3; 
H=k*(x*NH)^2/(2*9.8*m), J=0.5*k*(x*NH)^2

'pzt'
x=3.671e-6; f=1e-6; k = f/x, kpzt=k;
x=1.128e-5; J1 = 0.5 * k * x^2 / (200e-4 * 2e-4 * 4e-4);
Nm = k*x * x, L=200e-6; w=4e-6; h=2e-6; NW=cm/L; NH=cm/(w+2e-6); ND=cm/h; k=NW*k; k=k/NH; k=ND*k; m=7500*cm^3; 
H=k*(x*NH)^2/(2*9.8*m), J=0.5*k*(x*NH)^2

%'axial'
%x=3.543e-7; f=1000e-6; k=f/x;
%x=200e-6 * 0.15/100; J1 = 0.5 * k * x^2 / (200e-4 * 2e-4 * 4e-4);
%Nm = k*x * x, L=200e-6; w=4e-6; h=2e-6; NW=cm/(w+2e-6); NH=cm/L; ND=cm/h; k=NW*k; k=k/NH; k=ND*k; m=7500*cm^3; 
%H=k*(x*NH)^2/(2*9.8*m), J=0.5*k*(x*NH)^2
k
kpzt
g=9.8;h=0.5;k=2*m*g*h/(400e-6)^2
kk=kpzt*NW*NH*20

Nm = kk * (400e-6)^2
J = 0.5 * kk * (400e-6)^2 / (20*2e-4 * cm * cm)

%m=2300*cm^3 * 0.25 
%7500*cm*cm*5*2e-6

m=0.6e-3;g=9.8;h=1;k=2*m*g*h/(400e-6)^2

n = 270000 / 4 / 4 / (cm/200e-6) 
d = n * 6e-6

sheetlength = 270000 / (0.8*0.01 / 6e-6) * (202e-6)
%4 / 4 / (cm/200e-6) 
%d = n * 6e-6

wd=73.5e3 * (400e-6)^2 / 2 / (0.2*0.1 * cm * cm + 4*cm* 0.8*cm * 2e-4)

NM = 73.5e3 * (400e-6)^2

m = 73.5e3*10e-6 / 9.8

v = sqrt(73.5e3*(40e-6)^2 / 0.001)

Nm = 73.5e3*(400e-6)^2
