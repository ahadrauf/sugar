clear all;
E=160e9;
h=2e-6;
w=2e-6;
L=100e-6;
L3=20e-6;
I=w^3*h/12;
i=0;

T= 3.640;
Fy=10e-6;
mesh=100e-6;
Displacement_yt=12.1205e-6;
Displacement_yb=12.1203e-6;
Displacement_xt=-0.235603e-6; 
Displacement_xb=4.94275e-6;
i=i+1;
Ti(i)=T;
Fyi(i)=Fy;
meshi(i)=mesh;
L2=L3+sqrt((Displacement_yt-Displacement_yb)^2 + (Displacement_xt-Displacement_xb)^2);
alpha(i)=acos((Displacement_xt-Displacement_xb)/L2);
Displacement_yti(i)=Displacement_yt;
Displacement_xti(i)=Displacement_xt;
Displacement_ybi(i)=Displacement_yb;
Displacement_xbi(i)=Displacement_xb;
analytical_y(i)=Fy*L^3/3/E/I;
analytical_theta(i)=Fy*L^2/2/E/I;

T= 3.760;
Fy=20e-6;
mesh=100e-6;
Displacement_yt=24.2411e-6;
Displacement_yb=24.2406e-6;
Displacement_xt=-0.471205e-6; 
Displacement_xb=9.88549e-6;
i=i+1;
Ti(i)=T;
Fyi(i)=Fy;
meshi(i)=mesh;
L2=L3+sqrt((Displacement_yt-Displacement_yb)^2 + (Displacement_xt-Displacement_xb)^2);
alpha(i)=acos((Displacement_xt-Displacement_xb)/L2);
Displacement_yti(i)=Displacement_yt;
Displacement_xti(i)=Displacement_xt;
Displacement_ybi(i)=Displacement_yb;
Displacement_xbi(i)=Displacement_xb;
analytical_y(i)=Fy*L^3/3/E/I;
analytical_theta(i)=Fy*L^2/2/E/I;

T= 3.67;
Fy=30e-6;
mesh=100e-6;
Displacement_yt=36.3616e-6;
Displacement_yb=36.3609e-6;
Displacement_xt=-0.706808e-6;  
Displacement_xb=14.8282e-6;
i=i+1;
Ti(i)=T;
Fyi(i)=Fy;
meshi(i)=mesh;
L2=L3+sqrt((Displacement_yt-Displacement_yb)^2 + (Displacement_xt-Displacement_xb)^2);
alpha(i)=acos((Displacement_xt-Displacement_xb)/L2);
Displacement_yti(i)=Displacement_yt;
Displacement_xti(i)=Displacement_xt;
Displacement_ybi(i)=Displacement_yb;
Displacement_xbi(i)=Displacement_xb;
analytical_y(i)=Fy*L^3/3/E/I;
analytical_theta(i)=Fy*L^2/2/E/I;

T= 4.010;
Fy=40e-6;
mesh=100e-6;
Displacement_yt=48.4821e-6;
Displacement_yb=48.4812e-6;
Displacement_xt=-0.942411e-6; 
Displacement_xb=19.771e-6;
i=i+1;
Ti(i)=T;
Fyi(i)=Fy;
meshi(i)=mesh;
L2=L3+sqrt((Displacement_yt-Displacement_yb)^2 + (Displacement_xt-Displacement_xb)^2);
alpha(i)=acos((Displacement_xt-Displacement_xb)/L2);
Displacement_yti(i)=Displacement_yt;
Displacement_xti(i)=Displacement_xt;
Displacement_ybi(i)=Displacement_yb;
Displacement_xbi(i)=Displacement_xb;
analytical_y(i)=Fy*L^3/3/E/I;
analytical_theta(i)=Fy*L^2/2/E/I;

T= 4.040;
Fy=50e-6;
mesh=100e-6;
Displacement_yt=60.6026e-6;
Displacement_yb=60.6015e-6;
Displacement_xt=-1.17801e-6;  
Displacement_xb=24.7137e-6;  
i=i+1;
Ti(i)=T;
Fyi(i)=Fy;
meshi(i)=mesh;
L2=L3+sqrt((Displacement_yt-Displacement_yb)^2 + (Displacement_xt-Displacement_xb)^2);
alpha(i)=acos((Displacement_xt-Displacement_xb)/L2);
Displacement_yti(i)=Displacement_yt;
Displacement_xti(i)=Displacement_xt;
Displacement_ybi(i)=Displacement_yb;
Displacement_xbi(i)=Displacement_xb;
analytical_y(i)=Fy*L^3/3/E/I;
analytical_theta(i)=Fy*L^2/2/E/I;

%plots =============================================================================================

%Fy vs y
figure(1); clf; hold on; grid on;
plot(Fyi,Displacement_yti,'-ob');
plot(Fyi,analytical_y,'o-r');
title('Fy vs y, Linear beam, 2*2*100um^3, mesh=100um');
ylabel('y [N]   E=160GPa');
xlabel('Fy [m]    blue=IntelliSuite, Red=Roark');

%Fy vs time
figure(3); clf; hold on; grid on;
plot(Fyi,Ti,'-ob');
title('Fy vs t, Linear beam, 2*2*100um^3, mesh=100um');
xlabel('Fy [N]   blue=IntelliSuite');
ylabel('t [seconds]    E=160GPa');

%Fy vs relerr y
figure(5); clf; hold on; grid on;
Displacement_yti
analytical_y
plot(Fyi,(analytical_y-Displacement_yti)./analytical_y,'-og');
title('Fy vs RelErr(y), Linear beam, 2*2*100um^3, mesh=100um');
ylabel('RelativeError(y) [N]   E=160GPa');
xlabel('Fy [m]    blue=IntelliSuite, Red=Roark');

%Fy vs theta
figure(2); clf; hold on; grid on;
plot(Fyi(i:-1:2),-pi/2+alpha(i:-1:2),'-ob');
plot(Fyi(i:-1:2),analytical_theta(i:-1:2),'o-r');
title('Fy vs angle, Linear beam, 2*2*100um^3, mesh=100um');
ylabel('angle [rad]   E=160GPa');
xlabel('Fy [N]    blue=IntelliSuite, Red=Roark');
%Fy vs relerr theta
figure(4); clf; hold on; grid on;
theta2=analytical_theta(i:-1:2);
alpha2=-pi/2+alpha(i:-1:2);
plot(Fyi(i:-1:2),(theta2-alpha2)./theta2,'o-g');
title('Fy vs RelErr(angle), Linear beam, 2*2*100um^3, mesh=100um');
ylabel('RelativeError(angle) [rad]   E=160GPa');
xlabel('Fy [N]    blue=IntelliSuite, Red=Roark');

