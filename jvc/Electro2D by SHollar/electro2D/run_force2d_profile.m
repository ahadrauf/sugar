%run_force2d_profile
%Dimensions in microns, force is in microNewtons

clear all;
tic;
file1 = 'seg_a_gap_2.m';
file2 = 'seg_a_slot_1.m';   %force caculated on this one


seg1 = load_and_parse(file1);
seg2 = load_and_parse(file2);


%dimensions are in microns
G = 0.05;
L = 10;
H = 10;
W = 0.05;
Y=0;
z_depth = L ; % microns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii=0;jj=1;
%for Y = 2:-(2-0.1)/20:0.1
%for G = 2:-(2-0.1)/20:0.1
x11=0;      y11=0; x12=W;       y12=0;  x13=W;      y13=H;  x14=0;      y14=H;
x21=0+W+G;  y21=Y; x22=W+W+G;   y22=Y;  x23=W+W+G;  y23=Y+H;  x24=0+W+G;  y24=Y+H;
V1=50;
V2=-50;
seg1 = [x11 y11 x12 y12 V1
        x12 y12 x13 y13 V1
        x13 y13 x14 y14 V1
        x14 y14 x11 y11 V1];

seg2 = [x21 y21 x22 y22 V2
        x22 y22 x23 y23 V2
        x23 y23 x24 y24 V2
        x24 y24 x21 y21 V2];

delta_approx = 0.05 ; %microns



figure(1);plot_electro2d(seg1, seg2, delta_approx);grid on;

permitivity = 8.854e-6;% mks units * 1e-6 (avoids truncation error)

[F_xx, F_yy, tot_charge, mseg1, mseg2] = force2d_profile2(seg1,seg2,...
	 delta_approx, permitivity, z_depth);
 
size_seg1 = size(mseg1,1);
size_seg2 = size(mseg2,1);
 
% charge_1 is the charge vector for seg1
charge_1=[ones(size_seg1,1) ; zeros(size_seg2,1)].*tot_charge;
% charge_2 is the charge vector for seg2
charge_2=[zeros(size_seg1,1) ; ones(size_seg2,1)].*tot_charge;

% calculate total force (sum of force contributions from segment1
% onto segment2
% value is in micro newtons

%net force on segment 2
FF_X=charge_2'*F_xx'*charge_1;
FF_Y=charge_2'*F_yy'*charge_1;
 
j=0;
for i=1:length(charge_2)
    if charge_2(i)~=0
        j=j+1;
        q(j)=charge_2(i);
    end
end
figure(2); plot(1:j,q);grid on;
toc;

j=0;
for i=1:length(charge_1)
    if charge_1(i)~=0
        j=j+1;
        q(j)=charge_1(i);
    end
end
%figure(3); plot(1:j,q);grid on;

%%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
Fx(ii,jj)=FF_X;
Fy(ii,jj)=FF_Y;
gap(ii,jj)=G;
why(jj,jj)=Y;
%end
jj=jj+1
ii=0;
%end
%figure(4);plot(gap(:,1),abs(Fx(:,1)));grid on;


charge = sum(charge_2)
chargevssigma = sum(charge_2)*(H*1e-6*L*1e-6)
capacitanceQV = charge/V1
capacitanceEAD = 8.854e-12 * (H*1e-6)^2 / (G*1e-6)
Fcqv = charge^2 / (8.854e-12 * 2 * (H*1e-6*L*1e-6))
Fcead = capacitanceEAD^2 * V1^2 / (8.854e-12 * 2 * (H*1e-6*L*1e-6))
Fpp= 1/2 * 8.854e-12 * (H*1e-6*L*1e-6) * V1^2 / (G*1e-6)^2
G
H
L
FF_X
C = sqrt(abs(FF_X) / (V1^2 / (8.854e-12 * 2 * (H*1e-6*L*1e-6))))
Q=C*V1
