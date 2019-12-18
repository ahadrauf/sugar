%force.m
%Calculates the force on s2 due to all the s1's.

function force(s1,s2)

seg1 = load_and_parse(s1);
seg2 = load_and_parse(s2);

delta_approx = 0.05; 		%number of microns that the segments get divided.
permitivity = 8.854e-6;	%mks units * 1e-6 (avoids truncation error).
z_depth = 200; 				%microns, thickness.

[F_xx, F_yy, tot_charge, mseg1, mseg2] = force2d_profile(seg1,seg2,...
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
FF_X=charge_2'*F_xx'*charge_1
FF_Y=charge_2'*F_yy'*charge_1	
 
plot_electro2d(seg1,seg2, delta_approx);
