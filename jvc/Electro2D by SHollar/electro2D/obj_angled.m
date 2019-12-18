%figure(3);plot(phi*180/pi,F_X,'b',phi*180/pi,F_Y,'r');grid on;xlabel('angle from vertical [deg]');ylabel('blue = force_x, red = force_y [N]');title('Force vs angle for 2 plates. V=15V, xyz=2u*20u*20u');
%range=1:length(Q0{1});figure(3);plot(range,-Q1{1}(range),'b',range,-Q5{1}(range),'r',range,-Q0{1}(range),'k');grid on;xlabel('discretization [node]');ylabel('charge [C]');title('charge distribution at 0, 5, & 10 degrees. V=15V, xyz=2u*20u*20u');
%ds=0.1e-6;L=20e-6;F=1/(2*8.854e-12 * ds * L); range=1:length(Q0{1});lhs0;figure(3);plot(range,Q1{1}(range).^2 * F,'b',range,Q5{1}(range).^2 * F,'r',range,Q0{1}(range).^2 * F,'k');grid on;xlabel('discretization [node]');ylabel('force [N]');title('force distribution at 0, 5, & 10 degrees. V=15V, xyz=2u*20u*20u');
%V=-15;ds=0.1e-6;L=20e-6;F=1/(2*8.854e-12 * ds * L); range=1:length(Q0{1});bhs0;figure(3);plot(range,8.854e-12 * ds*L / 2e-6 * range./range,'g', range,Q1{1}(range)./V,'b',range,Q5{1}(range)./V,'r',range,Q0{1}(range)./V,'k');grid on;xlabel('discretization [node]');ylabel('capacitance [F]');title('capacitance distribution at 0, 5, & 10 degrees. V=15V, xyz=2u*20u*20u');
% sum(Q1{1}(range)./V) / (8.854e-12 * L^2 / 2e-6)
%figure(2); plot([rhs,221], 8.854e-12 * ds * L * V ./ Q1{1}([rhs,221]));
%rr=0; figure(2); r=21+rr:220-rr; plot([r], 8.854e-12 * ds * L * V ./ Q1{1}(r));percent=(1-rr*2/length(rhs))*100; hold on; grid on;
%[phi1,gap,F_X,F_Y,cap,total_charge_on_seg2,F_pp,capacitance_pp,Q_pp,Q1,rhs1,lhs1,bhs1,ths1]=obj_angled;
%range=1:length(Q5{1});figure(3);plot(range, -V * 8.854e-12 *ds*L / 2e-6 * range./range,'g', range,-Q1{1}(range),'b');grid on;xlabel('discretization [node]');ylabel('charge [C]');title('charge distribution at 0, 5, & 10 degrees. V=15V, xyz=2u*20u*20u');
%range=1:length(Q5{1});figure(3);plot(range,-Q0{1}(range),'b',range, -V * 8.854e-12 *ds*L / 2e-6,'g');grid on;xlabel('discretization [node]');ylabel('charge [C]');title('charge distribution at 0, 5, & 10 degrees. V=15V, xyz=2u*20u*20u');
%FL=sum( Q0{1}(ceil([lhs]))  )^2 /2/8.854e-12/(20e-6 * 20e-6)
%FR=sum( Q0{1}(ceil([rhs]))  )^2 /2/8.854e-12/(20e-6 * 20e-6)
%PPSC_ratio=Qppsc/sum(Q{1}(rhs)),PP_ratio=Qpp/sum(Q{1}(rhs)), dQppsc,toc

function [phi,gap,F_X,F_Y,cap,totQ,F_pp,capacitance_pp,Qpp,Q,rhs,lhs,bhs,ths,Qppsc,dQppsc,deltaQ]=obj_angled
%init
D = [];

%dimensions [microns]
%Voltage [volts]
w1=1; 
h1= 20;
L = 20;
V1=15;V2=0;
ds = 0.1 + eps;

j = 0;
%for G = 2: 8/20  : 10
%for ang = 0 : 0.05 : 1
G=1; %gap
D = [];
j = j + 1;
x20 = w1 + G;
gap(j) = G;
ang = 0;
theta = ang * pi/18;
phi(j) = theta;

%Define rectangles here.
%The first one will be assigned to seg2, which is the one that is analized
%[x,y,w,h,V]

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
a{1} = [   0,    0,      w1,     h1,     V1]; %seg2. The one analyzed. 
a{2} = [ x20,    0,      w1,     h1,     V2];
A = [a{1}; a{2}];

for i = 1:size(A,1)
    x0=a{i}(1);
    y0=a{i}(2);
    w=a{i}(3);
    h=a{i}(4);
    V=a{i}(5);
    xa=x0-w/2; ya=y0-h/2;
    xb=x0+w/2; yb=y0-h/2;
    xc=x0+w/2; yc=y0+h/2;
    xd=x0-w/2; yd=y0+h/2;
    if i==1
        rot = R * [xa xb xc xd; ya yb yc yd];
        xa = rot(1,1); ya = rot(2,1);
        xb = rot(1,2); yb = rot(2,2);
        xc = rot(1,3); yc = rot(2,3);
        xd = rot(1,4); yd = rot(2,4);
    end
    C = [xa ya xb yb V; xb yb xc yc V; xc yc xd yd V; xd yd xa ya V]; 
    D = [D;C];
end

seg2 = D(1:4,:);
seg1 = D(5:size(D,1),:);

[charge_1, charge_2, capacitance, FF_X, FF_Y] = electro2d_2(seg1, seg2, ds, L);

cap(j) = capacitance;
F_X(j) = FF_X*1e-6;
F_Y(j) = FF_Y*1e-6;
totQ(j) = 1e-12 * L * ds * abs(sum(charge_2));

H=h1; %G=0.1;
F_pp(j) = 1/2 * 8.854e-12 * (H*1e-6 * L*1e-6) * V1^2 / (G*1e-6)^2;
capacitance_pp(j) = 8.854e-12 * (H*1e-6)*(L*1e-6) / (G*1e-6);
Qpp(j) = capacitance_pp(j)*abs(V1-V2);

w=w1;h=h1;
bhs = 1:w/ds;
rhs = w/ds+1:(w+h)/ds;
ths = (w+h)/ds+1:(2*w+h)/ds;
lhs = (2*w+h)/ds+1:(2*w+2*h)/ds;
Q{j} = 1e-12 * L * ds * charge_2(find(charge_2 ~= 0));

%end

figure(1);plot_electro2d(seg1, seg2, ds);grid on; axis equal; %objects
%figure(2); plot(1:length(rhs) ,Q(rhs));grid on; %charge on seg2
%figure(3); plot(1:length(Q),Q);grid on; %charge on rhs

dQppsc = Qpp/(L*1e-6 * h1*1e-6) * (G*1e-6 / pi) * L*1e-6;
Qppsc = dQppsc*2 + Qpp;

deltaQ=sum( -Q{1}([ceil(rhs),ceil(max(rhs)+1)]) - Qpp/(length(rhs)+1) );
