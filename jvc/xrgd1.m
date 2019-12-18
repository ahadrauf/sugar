%first beam
L1=100e-6;
b1i=-L1*sin(pi/4);
a1i=-L1*cos(pi/4);
b1j=0;
a1j=L1;
tau1i = [1 0 0; 0 1 0; b1i -a1i 1];
tau1j = [1 0 0; 0 1 0; b1j -a1j 1];
o3=zeros(3);
T10=[tau1i o3;o3 tau1j];
kij = [...
  6.6000e+003            0            0 -6.6000e+003            0            0
            0  2.6400e+000  1.3200e-004            0 -2.6400e+000  1.3200e-004
            0  1.3200e-004  8.8000e-009            0 -1.3200e-004  4.4000e-009
 -6.6000e+003            0            0  6.6000e+003            0            0
            0 -2.6400e+000 -1.3200e-004            0  2.6400e+000 -1.3200e-004
            0  1.3200e-004  4.4000e-009            0 -1.3200e-004  8.8000e-009 ];
T1=rot2local(0,0,pi/2)';
T1=[T1 o3;o3 T1];
k1ij=T1*kij*T1';
k1qr=T10*k1ij*T10';

%second beam
L2=100e-6;
b2i=0;
a2i=-L2;
b2j=L2;
a2j=0;
tau2i = [1 0 0; 0 1 0; b2i -a2i 1];
tau2j = [1 0 0; 0 1 0; b2j -a2j 1];
o3=zeros(3);
T20=[tau2i o3;o3 tau2j];
kij = [...
  6.6000e+003            0            0 -6.6000e+003            0            0
            0  2.6400e+000  1.3200e-004            0 -2.6400e+000  1.3200e-004
            0  1.3200e-004  8.8000e-009            0 -1.3200e-004  4.4000e-009
 -6.6000e+003            0            0  6.6000e+003            0            0
            0 -2.6400e+000 -1.3200e-004            0  2.6400e+000 -1.3200e-004
            0  1.3200e-004  4.4000e-009            0 -1.3200e-004  8.8000e-009 ];
T2=rot2local(0,0,pi)';
T2=[T2 o3;o3 T2];
k2ij=T2*kij*T2';
k2qr=T20*k1ij*T20';

%assemble
K1=zeros(9,9);
K2=zeros(9,9);
K1(1:6,1:6)=k1qr;
K2(4:9,4:9)=k2qr;
KQR=K1+K2;

%remove anchor
KQR(:,1)=[];KQR(:,1)=[];KQR(:,1)=[];
KQR(1,:)=[];KQR(1,:)=[];KQR(1,:)=[];

F=0*5e-6;
UQR=inv(KQR)*[0 0 0, F 0 0]';
UQR=[0;0;0;UQR];

UIJ=[T10'*UQR(1:6,1) ; T20'*UQR(4:9,1)];

u=[UQR(1:3,1);UIJ(1:3,1);UIJ(4:6,1);UQR(4:6,1);UIJ(7:9,1);UIJ(10:12,1);UQR(7:9,1)];
pos{1}=[0;0;0];
x=L1*cos(pi/4);
y=L1*sin(pi/4);
pos{2}=[x;y;0];
x=x+L1*cos(-pi/4);
y=y+L1*sin(-pi/4);
pos{3}=[x;y;0];
x=x+L1;
y=y+0;
pos{4}=[x;y;0];
x=x+L2;
y=y+0;
pos{5}=[x;y;0];
x=x+0;
y=y+L2;
pos{6}=[x;y;0];

R{1}=rot2local(0,0,pi/4);
R{2}=rot2local(0,0,-pi/4);
R{3}=rot2local(0,0,0);
R{4}=rot2local(0,0,0);
R{5}=rot2local(0,0,pi/2);
R{6}=rot2local(0,0,pi/2);
l=L1;
w(1)=8e-6;
w(2)=2e-6;
w(3)=8e-6;
w(4)=8e-6;
w(5)=2e-6;
w(6)=8e-6;
h=2e-6;

figure(1);
clf;
grid on;
hold on;
view(0,90);
xlabel('X - horizontal  [m]');		%x-axis label.
ylabel('Y - vertical  [m]');		%y-axis label.
zlabel('Z - out of plane  [m]');	%z-axis label.

   Q=zeros(12,1);
   displaybeam(Q, pos{1}, rot2local(0,0,pi)', 20e-6, 20e-6, h);
for i=1:6 
   Q=zeros(12,1);
   q=u(i*3 - 2:i*3 -2 + 5,1);
   Q([1 2 6 7 8 12]) = q(1:6,1);
   displaybeam(Q, pos{i}, R{i}', l, w(i), h);
end

rotate3d on;
color=pink;
color(:,2)=color(:,2)*(0.5 + 0.5*rand);
color(:,1)=color(:,1)*(0.5 + 0.5*rand);
color(:,3)=color(:,3)*(0.5 + 0.5*rand);
colormap(color);
shading interp;
axis equal;
axis vis3d;
hold off;
