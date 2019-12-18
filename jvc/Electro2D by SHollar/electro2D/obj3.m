%[gap,F_X,F_Y,cap,total_charge_on_seg2,F_pp,capacitance_pp,Q_pp,Q,rhs,lhs,bhs,ths]=obj3


function [OL,F_X,F_Y,cap,total_charge_on_seg2,F_pp,cap_pp,Q_pp,Q,rhs,lhs,bhs,ths]=obj_slot
%init
D = [];

%dimensions [microns]
%Voltage [volts]
w1=4; h1=2;
w2=w1;
V1=50;V2=0;
L = 100;
ds = 0.5 + eps;

j = 0;
for overlap = 0: 0.25 : 5
%overlap=10
D = [];
j = j + 1;
x20 = overlap;
OL(j) = overlap;

%Define rectangles here.
%The first one will be assigned to seg2, which is the one that is analized
%[x,y,w,h,V]
%g=1.375;
g=0.5;
gg=0.1;
a{1} = [overlap,      4,      4,     h1,     V2]; %seg2. The one analyzed. 
a{2} = [0,            4-2-g,      w1,     h1,     V2]; %lower
a{3} = [0,            4+2+g,      w1,     h1,     V2];  %upper
a{4} = [0,            4-2-g-2-gg,      w2,     h1,     V1];  %lower-lower
a{5} = [0,            4+2+g+2+gg,      w2,     h1,     V1];  %upper-upper
A = [a{1}; a{2}; a{3}; a{4}; a{5}];

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
    C = [xa ya xb yb V; xb yb xc yc V; xc yc xd yd V; xd yd xa ya V];
    D = [D;C];
end

seg2 = D(1:4,:);
seg1 = D(5:size(D,1),:);

[charge_1, charge_2, capacitance, FF_X, FF_Y] = electro2d_2(seg1, seg2, ds, L);

cap(j) = capacitance;
F_X(j) = FF_X*1e-6
F_Y(j) = FF_Y*1e-6
total_charge_on_seg2(j) = 1e-12 * L * ds * abs(sum(charge_2));

G=2;
F_pp(j) = 2 * 1/2 * 8.854e-12 * (L*1e-6) * V1^2 / (G*1e-6);
if overlap > 10
    over=0;
else
    over=overlap;
end
cap_pp(j) = 2 * 8.854e-12 * (over*1e-6)*(L*1e-6) / (G*1e-6);
Q_pp(j) = cap_pp(j)/abs(V1-V2);

w=w1;h=h1;
bhs = 1:w/ds;
rhs = w/ds+1:(w+h)/ds;
ths = (w+h)/ds+1:(2*w+h)/ds;
lhs = (2*w+h)/ds+1:(2*w+2*h)/ds;
Q{j} = charge_2(find(charge_2 ~= 0));

end

figure(1);plot_electro2d(seg1*1e-6, seg2*1e-6, ds*1e-6);grid off; axis equal; %objects
%figure(2); plot(1:length(rhs) ,Q(rhs));grid on; %charge on seg2
%figure(3); plot(1:length(Q),Q);grid on; %charge on rhs
