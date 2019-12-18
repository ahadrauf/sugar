function [capacitance, FF_X, FF_Y] = electro2d(seg1, seg2, ...
	 delta_approx, permitivity, z_depth);

%function [capacitance, FF_X, FF_Y] = electro2d(seg1, seg2, ...
%	   delta_approx, permitivity, z_depth);
%
%        electro2d is function that models the 2d electrstatic 
%        behaviour of multiple conducting bodies.  The model is 
%        based on the theory of a quasi-static electric field 
%        derived from Maxwell's equations.  Assuming the 
%        surfaces of the conducting bodies are specified in seg1 
%        and seg2, electro2d calculates the capacitance and forces 
%        acting on seg2 with respect to seg1.
%        
%        capacitance (Farads) is the capacitance between the 
%        bodies as defined by seg1 and seg2.  capacitance is ONLY
%        valid when all segments in seg1 have the same voltage and
%        all segments in seg2 have the same voltage.
%
%        FF_X (microNewtons) = force acting on seg2 in the x direction
%        due to seg1
% 
%        FF_Y (microNewtons) = force acting on seg2 in the y direction
%        due to seg1
%
%        seg1 (microns, Volts) = array of line segments defining 
%        a conducting body.
%        Each line is defined as [x1 y1 x2 y2 Voltage] where  (x1,y1) 
%        and (x2,y2) define the endpoints of the line and Voltage defines
%        the voltage on the line.
%
%        seg2 = same as seg1
%
%        delta_approx (microns) = seg1 and seg2 are broken up 
%        into smaller line segments to approximate the ideal 
%        behavior.  The length of these smaller line segments 
%        is delta approx.  The smaller the delta_approx the
%        more accurate the results.  Nevertheless, simulation 
%        time is on order O(N^2) where N=1/delta_approx (ie 
%        the smaller the delta approx, the longer it takes 
%        to simulate).
%
%        permitivity of free space is 8.854e-6 (mks units) * 1e-6.  
%        Permitivity value is scaled by 1e6 to avoid truncation error
%        in matlab.
%
%        z_depth = since this is only a 2D model, z_depth defines the third
%        dimension (ie the depth of conducting bodies out of the plane) 
%
%        see also plot_electro2d

epsilon = permitivity; % program calculates based on 1e6 scale factor
delta=1e-15; % truncation error factor
z_axis=z_depth;%microns


% this section takes in seg1 and seg2 and breaks the lines up 
% into smaller line segments of length delta approx

adding_matrix=[0 0 0 0];
for	n=1:size(seg1,1)
   adding_matrix=[adding_matrix;dissect(seg1(n,:),delta_approx)];
end
volt_m1=adding_matrix(2:size(adding_matrix,1),:);

adding_matrix=[0 0 0 0];
for	n=1:size(seg2,1)
   adding_matrix=[adding_matrix;dissect(seg2(n,:),delta_approx)];
end
volt_m2=adding_matrix(2:size(adding_matrix,1),:);

size_volt2 = size(volt_m2,1);
size_volt1 = size(volt_m1,1);

% total segment number represents the total number of small segments

total_segment_number=size_volt1+size_volt2


'matrix calculations'

% all calculations are done in matrix form.  This consumes quite a bit
% of memory but since matlab does fast matrix calculations, it is 
% very fast

% rotate line segments so that they lie on the x axis.  This is needed
% to do the force calculations

del_theta=ones(size_volt1,1)*volt_m2(:,3)'...
   -volt_m1(:,3)*ones(1,size_volt2);
costh=cos(del_theta);
sinth=sin(del_theta);
clear del_theta

rot_theta=volt_m1(:,3)*ones(1,size_volt2);
cos_rot_theta=cos(rot_theta);
sin_rot_theta=sin(rot_theta);
clear rot_theta

del_x=ones(size_volt1,1)*volt_m2(:,1)'...
   -volt_m1(:,1)*ones(1,size_volt2);
del_y=ones(size_volt1,1)*volt_m2(:,2)'...
   -volt_m1(:,2)*ones(1,size_volt2);

X_trans=cos_rot_theta.*del_x...
   +sin_rot_theta.*del_y;
Y_trans=cos_rot_theta.*del_y...
   -sin_rot_theta.*del_x;
clear del_x del_y

Y_cos=Y_trans.*costh;
R_par=X_trans.*costh+Y_trans.*sinth;
R_perp=-X_trans.*sinth+Y_cos;
div_sinth=(2*(sinth>=0)-1)./(abs(sinth)+delta);


% once rotated we can now calculate the force contribution from each segment
% in seg1 to each segment in seg2.  Matrices are of the size 
% size_volt_m1 X size_volt_m2

'calculate force matrices'

L=delta_approx/2;
denom_plusL=Y_trans+L*sinth;
denom_minusL=Y_trans-L*sinth;
less_zero_plusL=(2*(denom_plusL>=0)-1);
less_zero_minusL=(2*(denom_minusL>=0)-1);

alpha=-delta_approx/2;
sub_force_matrices;
sub_x_sum=sub_x;
sub_y_sum=sub_y;
y_force_discon=-y_discont;

alpha= delta_approx/2;

% sub_force_matrices in a .m file that actually calculates the force
% on each segment

sub_force_matrices;
sub_x_sum=sub_x_sum-sub_x;
sub_y_sum=sub_y_sum-sub_y;
y_force_discon=y_force_discon+y_discont;

sub_y_sum=sub_y_sum+y_force_discon...
   .*((less_zero_plusL.*less_zero_minusL)<0)...
   .*less_zero_plusL;

% include epsilon plus factors here

sub_x_sum=-1/2/pi/epsilon*sub_x_sum;
sub_y_sum=-1/2/pi/epsilon*sub_y_sum;

% now we have to rotate back to original position

F_xx=cos_rot_theta.*sub_x_sum-sin_rot_theta...
   .*sub_y_sum;
F_yy=sin_rot_theta.*sub_x_sum+cos_rot_theta...
   .*sub_y_sum;

% multiple by depth to scale in 3rd dimension
% F_xx is a force matrix from every line segment in seg1
% to every line segment in seg2.  When F_xx is multiplied by
% charge densities from seg1 and seg2, the total force is obtained.
% ie. total force in x = (charge density for seg1)*F_xx*(charge density
% for seg2)

F_xx=F_xx*z_axis;
F_yy=F_yy*z_axis;

% clear up some memory

clear cos_rot_theta sin_rot_theta sub_x_sum ...
  X_trans Y_trans sinth costh ...
  div_sinth R_par R_perp denom_plusL ...
  denom_minusL less_zero_plusL less_zero_minusL ...
  sinth costh Y_cos

'calculate potential matrix'

% now we have to determine the charge density on each line segment.
% we now the potential on each line segment (initial condition).
% What we do now is similar to the above force calculation.  We find
% the potential contribution from each line segment to every other line
% segment.  The matrices here are even larger (total_line_segments X 
% total_line_segments).

volt_m = [volt_m1 ; volt_m2];	
size_volt=size(volt_m,1);

% again we have to rotate the frame such that the line segment is 
% parallel to the X axis and centered at the origin

x_value=volt_m(:,1)*ones(1,size_volt)...
   -ones(size_volt,1)*volt_m(:,1)';
y_value=volt_m(:,2)*ones(1,size_volt)...
   -ones(size_volt,1)*volt_m(:,2)';
pot_angle=volt_m(:,3)*ones(1,size_volt);
sin_angle=sin(pot_angle);
cos_angle=cos(pot_angle);
x_pot=x_value.*cos_angle+y_value.*sin_angle;
y_pot=y_value.*cos_angle-x_value.*sin_angle;

x_p=x_pot-delta_approx/2;
x_m=x_pot+delta_approx/2;

% below is the actual calculation of the potnential contribution
% from one segment to another.

s_matrix=(x_m/2.*log(x_m.^2+y_pot.^2+delta)...
   		-x_p/2.*log(x_p.^2+y_pot.^2+delta) ...
   +y_pot.*atan(x_m./(y_pot+delta))...
   -y_pot.*atan(x_p./(y_pot+delta))...
   +delta_approx)/2/pi/epsilon;

s_matrix=s_matrix';

% s_plus_matrix adds an extra line.  This essentiall forces the sum
% of the charge of the entire system to be zero.

s_plus_matrix=[s_matrix ones(size_volt,1)];
s_plus_matrix=[s_plus_matrix; ones(1,size_volt+1)];
s_plus_matrix(size_volt+1,size_volt+1)=0;
voltage=volt_m(:,4);
voltage(size_volt+1)=0;

'invert potential matrix'
% since voltage = s_plus_matrix*sigma_charge
% we can invert the matrix to find the charge density
% (sigma_charge)

sigma_charge=inv(s_plus_matrix)*voltage;
   
%volt_m=[volt_m sigma_charge(1:size_volt)];  

clear s_plus_matrix s_matrix cos_angle ...
	sin_angle pot_angle x_value y_value ...
	x_pot y_pot x_p x_m

% remove last line of sigma_charge (doesn't represent charge)
tot_charge=sigma_charge(1:size(sigma_charge,1)-1);
% charge_1 is the charge vector for seg1
charge_1=tot_charge(1:size_volt1);
% charge_2 is the charge vector for seg2
charge_2=tot_charge(size_volt1+1:size(tot_charge,1));
% calculate total force (sum of force contributions from each line
% segment
% value is in micro newtons
FF_X=charge_2'*F_xx'*charge_1;
FF_Y=charge_2'*F_yy'*charge_1;	
mag_F=sqrt(FF_X^2+FF_Y^2);
% capacitance is defined as charge on seg2 divided by difference in voltage
% NOTE this only works when seg1 is all the same voltage and 
% seg2 is all the same voltage 
capacitance=1e-12*z_axis*delta_approx...
   *abs(sum(charge_2)/(seg1(1,5)-seg2(1,5) + eps));
























