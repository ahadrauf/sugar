
%dimensions are in microns
G = 0.1; %gap
L = 10;  %z-axis
H = 10;  %y-axis
W = 0.2; %x-axis
Y = 0;   %y-offset
z_depth = L;

x11=0;      y11=0;      x12=W;       y12=0;         x13=W;      y13=H;          x14=0;      y14=H;
x21=0+W+G;  y21=Y;      x22=W+W+G;   y22=Y;         x23=W+W+G;  y23=Y+H;        x24=0+W+G;  y24=Y+H;

V1=50; V2=0;

seg1 = [x11 y11 x12 y12 V1
        x12 y12 x13 y13 V1
        x13 y13 x14 y14 V1
        x14 y14 x11 y11 V1];

seg2 = [x21 y21 x22 y22 V2
        x22 y22 x23 y23 V2
        x23 y23 x24 y24 V2
        x24 y24 x21 y21 V2];

delta_approx = 0.1 ; %microns



[charge_1,charge_2,capacitance, FF_X, FF_Y] = electro2d_2(seg1, seg2, delta_approx, z_depth);
figure(1);plot_electro2d(seg1, seg2, delta_approx);grid on;

capacitance
F_X = FF_X*1e-6
F_Y = FF_Y*1e-6
total_charge_on_seg2 = 1e-12 * z_depth * delta_approx * abs(sum(charge_2))

F_pp= 1/2 * 8.854e-12 * (H*1e-6 * L*1e-6) * V1^2 / (G*1e-6)^2
capacitance_pp = 8.854e-12 * (H*1e-6)*(L*1e-6) / (G*1e-6)
Q_pp=capacitance_ead*abs(V1-V2)

Q=charge_2(find(charge_2 ~= 0));
figure(2); plot((1:length(Q))*delta_approx ,Q);grid on;

%right = (1+W/delta_approx) : (1+(H+W)/delta_approx);
%bottom = 1 : W/delta_approx;
%top = (1+(H+W)/delta_approx) : (1+(H+W)/delta_approx) + (1+W/delta_approx);
%figure(3); plot( 1:H/delta_approx, Q(1+W/delta_approx:(1+(H+W)/delta_approx)));

