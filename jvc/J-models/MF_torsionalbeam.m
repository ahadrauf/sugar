%R is a 3x3 rotation to local
%q is a 24x1 vector of displacement node1, displacement node2, velocity node1, velocity node2.
function [output] = MF_torsionalbeam(flag, R, param, q, t, nodes, varargin);

switch(flag)

case 'vars'

  output.dynamic = {1 {'x' 'y' 'z' 'rx' 'ry' 'rz'};
                    2 {'x' 'y' 'z' 'rx' 'ry' 'rz'}};

case 'check'

  if (~isfield(param, 'density')       | ...
      ~isfield(param, 'fluid')         | ...
      ~isfield(param, 'viscosity')     | ...
      ~isfield(param, 'Youngsmodulus'))
    output = 'Missing params typically specified in process file';
  elseif ~isfield(param, 'l')
    output = 'Missing length';
  elseif ~isfield(param, 'w')
    output = 'Missing width';
  elseif ~isfield(param, 'h') 
    output = 'Missing height (layer thickness)';
  else
    output = [];
  end


case 'M'
   
   %geometrical params
   H = param.h;               %Beam layer thickness, along z-axis.
   W = param.w;               %Beam width, along y-axis.
   L = param.l;               %Beam length, along x-axis.
   D = param.density;         %Poly-Si density.
   A = H * W;                 %Beam crossectional area, yz-plane.
   Iy = W * H^3 / 12;         %Rotatory moment of inertia about y.
   Iz =	H * W^3 / 12;         %Rotatory moment of inertia about z. Used in 2D.
   Jx =	Iy + Iz;              %Polar moment of inertia, torsional inertia.
                              %J=0.14 if h<=w, cJ=1/3 if h<<w.
   A = H * W;					%Cross sectional area.   
   a =  1 / 3;					%Various quantities...
   b = 13 / 35;
   c = 13 / 35;
   d = Jx / (3 * A);
   e = L^2 / 105;
   f = L^2 / 105;
   g = 1 / 3;
   h = 13 / 35;
   i = 13 / 35;
   j = Jx / (3 * A);
   k = L^2 / 105;
   l = L^2 / 105;
   m = -11 * L / 210;
   n = 11 * L / 210;
   oo= 	9 / 70;
   p = 13 * L / 420;
   q = 9 / 70;
   r = -13 * L / 420;
   s = Jx / (6 * A);
   t = -13 * L / 420;
   u = -L^2 / 140;
   v = 1 / 6;
   w = -L^2 / 140;
   
%rotation matrix
   Rall(1:3,1:3) = R;
   Rall(4:6,4:6) = R;
   R = Rall;
  	 
%mass matrix 
  
   M11 = R * [ ...
     a     0     0     0     0     0;
     0     b     0     0     0     n;
     0     0     c     0     m     0;
     0     0     0     d     0     0;
     0     0     m     0     e     0;     
     0     n     0     0     0     f] * R';
   M22 = R * [ ...
     g     0     0     0     0     0;
     0     h     0     0     0    -n;
     0     0     i     0    -m     0;
     0     0     0     j     0     0;
     0     0    -m     0     k     0;
     0    -n     0     0     0     l] * R';
   M21 = R * [ ...
     v     0     0     0     0     0;
     0     oo    0     0     0     p;
     0     0     q     0     r     0;
     0     0     0     s     0     0;
     0     0    -r     0     w     0;
     0     t     0     0     0     u] * R';
  
   output = D*A*L*[M11  M21'; M21 M22];
  
case 'D'

%matrix params
   H = param.h;		%Beam layer thickness, along z-axis.
   W = param.w;		%Beam width, along y-axis.
   L = param.l;		%Beam length, along x-axis.
   A = H * W;		%Beam crossectional area, yz-plane.
   Mu = param.viscosity;%Viscocity of fluid environment.
   delta = param.fluid;	%Fluid layer thickness.
   Iy =	W * H^3 / 12;	%Rotatory moment of inertia about y.
   Iz =	H * W^3 / 12;	%Rotatory moment of inertia about z. Used in 2D.
   Jx =	Iy + Iz;	%Polar moment of inertia, torsional inertia.
                        %J=0.14 if h<=w, cJ=1/3 if h<<w.
   A = H * W;					%Cross sectional area.   
   a = 1 / 3;					%Various quantities...
   b = 13 / 35;
   c = 13 / 35;
   d = Jx / (3 * A);
   e = L^2 / 105;
   f = L^2 / 105;
   g = 1 / 3;
   h = 13 / 35;
   i = 13 / 35;
   j = Jx / (3 * A);
   k = L^2 / 105;
   l = L^2 / 105;
   m = -11 * L / 210;
   n = 11 * L / 210;
   oo = 9 / 70;
   p = 13 * L / 420;
   q = 9 / 70;
   r = -13 * L / 420;
   s = Jx / (6 * A);
   t = -13 * L / 420;
   u = -L^2 / 140;
   v = 1 / 6;
   w = -L^2 / 140;
   
%rotation matrix
   Rall(1:3,1:3) = R;
   Rall(4:6,4:6) = R;
   R = Rall;
   
%damping matrix. This viscous model is good for planar motion.  
   D11 = R * [ ...
     a     0     0     0     0     0;
     0     b     0     0     0     n;
     0     0     c     0     m     0;
     0     0     0     d     0     0;
     0     0     m     0     e     0;     
     0     n     0     0     0     f] * R';
   D22 = R * [ ...
     g     0     0     0     0     0;
     0     h     0     0     0    -n;
     0     0     i     0    -m     0;
     0     0     0     j     0     0;
     0     0    -m     0     k     0;
     0    -n     0     0     0     l] * R';
   D21 = R * [ ...
     v     0     0     0     0     0;
     0     oo    0     0     0     p;
     0     0     q     0     r     0;
     0     0     0     s     0     0;
     0     0    -r     0     w     0;
     0     t     0     0     0     u] * R';
  
  output = Mu*W*L/delta*[D11  D21'; D21 D22];

case 'K'

%matrix params
   E = param.Youngsmodulus;	%Modulus of elasticity.
   H = param.h;			%Beam layer thickness, along z-axis.
   W = param.w; 		%Beam width, along y-axis.
   L = param.l;			%Beam length, along x-axis.
   A = H * W;			%Beam crossectional area, yz-plane.

   Nu = param.Poisson;		%Poisson's ratio.
   G = E / (2 * (1 + Nu));	%Shear modulus. GJ = Torsional stiffness 

   if (W>H)
     J = W*H^3*(16/3-3.36*H/W*(1-(H/W)^4/12))/16;
   else
     J = H*W^3*(16/3-3.36*W/H*(1-(W/H)^4/12))/16;
   end
 
% The newer version used this different formula for J.  If the two
% are supposed to be the same, there is a very high relative error
% (at least ~0.9 for the mirror mode demo).  I don't know which one
% is correct.

%   J2 = 2/7*(W*H)^3/(W^2+H^2);	%Polar second moment of area.
%   abs((J2-J)/J)

   Iy =	W * H^3 / 12;		%Moment about y.
   Iz =	H * W^3 / 12;		%Moment about z. Used in 2D.

   a = E * A / L;
   b = 12 * E * Iz / L^3;
   c = 12 * E * Iy / L^3;
   d = G * J / L;
   e = 4 * E * Iy / L;
   f = 2 * E * Iy / L;
   g = 4 * E * Iz / L;
   h = 2 * E * Iz / L;
   i = 6 * E * Iy / L^2;
   j = 6 * E * Iz / L^2;
   
%rotation matrix
   Rall(1:3,1:3) = R;
   Rall(4:6,4:6) = R;
   R = Rall;
   
%stiffness matrix, linear model
   K11 = R * [ ...
     a     0     0     0     0     0; 
     0     b     0     0     0     j; 
     0     0     c     0    -i     0; 
     0     0     0     d     0     0; 
     0     0    -i     0     e     0;     
     0     j     0     0     0     g] * R';
   K22 = R * [ ...
     a     0     0     0     0     0; 
     0     b     0     0     0    -j; 
     0     0     c     0     i     0; 
     0     0     0     d     0     0; 
     0     0     i     0     e     0; 
     0    -j     0     0     0     g] * R';   
   K12 = R * [ ...         
    -a     0     0     0     0     0; 
     0    -b     0     0     0     j; 
     0     0    -c     0    -i     0; 
     0     0     0    -d     0     0; 
     0     0     i     0     f     0;     
     0     -j    0     0     0     h] * R';
  
   output = [K11  K12; K12' K22];


case 'pos'

  % Compute relative positions of beam nodes

  output = R * [0 param.l;
                0 0;
                0 0];            

case 'F'
%restoring force correction for torsional twist
   %Method: substract off the linear contribution, and add the nonlinear contribution 
   W=param.w; %width
   H=param.h; %layer thickness
   L=param.l; %beam length
   o3=zeros(3,3); %bunch of zeros
   o6=zeros(6,6); %bunch of zeros
   R6=[R o3; o3 R]; %6x6 rotation
   R12=[R6 o6; o6 R6]; %12x12 rotation
   E=param.Youngsmodulus; %Young's modulus.
   qloc=R12*q(1:12,1); %local displacements of node1
   y1=qloc(2); %y translation displacement node 1
   z1=qloc(3); %z translation displacement node 1
   ox1=qloc(4); %ox twist displacement node 1
   oz1=qloc(6); %z rotation displacement node 1
   oy1=qloc(5); %y rotation displacement node 1
   y2=qloc(8); %y translation displacement node 2
   z2=qloc(9); %z translation displacement node 2
   ox2=qloc(10); %ox twist displacement node 2
   oz2=qloc(12); %z rotation displacement node 2
   oy2=qloc(11); %y rotation displacement node 2
   Iz=H*W^3/12; %Moment of inertial about z-axis.
   Iy=H^3*W/12; %Moment of inertial about y-axis.
   E12LLL=E*12/L^3; %part of the stiffness, see K.
   E6LL=E*6/L^2; %part of the stiffness, see K.
   
   nonlinFy1=(y1-y2)*E12LLL*( Iz/2*(1+cos(ox1)) + Iy/2*sin(ox1) ) + (oz1+oz2)*E6LL*( Iz/2*(1+cos(ox1)) + Iy/2*sin(ox1) ); %y nonlinear stiffness including torsion node 1
   nonlinFy2=(y2-y1)*E12LLL*( Iz/2*(1+cos(ox2)) + Iy/2*sin(ox2) ) - (oz1+oz2)*E6LL*( Iz/2*(1+cos(ox2)) + Iy/2*sin(ox2) ); %y nonlinear stiffness including torsion node 2
   linFy1=(y1-y2)*E12LLL*Iz + (oz1+oz2)*E6LL*Iz; %y linear stiffness node 1
   linFy2=(y2-y1)*E12LLL*Iz - (oz1+oz2)*E6LL*Iz; %y linear stiffness node 2
   
   nonlinFz1=(z1-z2)*E12LLL*( Iy/2*(1+cos(ox1)) + Iz/2*sin(ox1) ) + (oz1+oz2)*E6LL*( Iy/2*(1+cos(ox1)) + Iz/2*sin(ox1) ); %z nonlinear stiffness including torsion node 1
   nonlinFz2=(z2-z1)*E12LLL*( Iy/2*(1+cos(ox2)) + Iz/2*sin(ox2) ) - (oz1+oz2)*E6LL*( Iy/2*(1+cos(ox2)) + Iz/2*sin(ox2) ); %z nonlinear stiffness including torsion node 2
   linFz1=(z1-z2)*E12LLL*Iy + (oy1+oy2)*E6LL*Iy; %z linear stiffness node 1
   linFz2=(z2-z1)*E12LLL*Iy - (oy1+oy2)*E6LL*Iy; %z linear stiffness node 2
   
   
   nonlinFy1=(y1-y2)*E12LLL*( Iz/2*(1-cos(ox1)) + Iy/2*(1-cos(ox1)) ) + (oz1+oz2)*E6LL*( Iz/2*(1-cos(ox1)) + Iy/2*(1-cos(ox1)) ); %y nonlinear stiffness including torsion node 1
   nonlinFy2=(y2-y1)*E12LLL*( Iz/2*(1-cos(ox2)) + Iy/2*(1-cos(ox2)) ) - (oz1+oz2)*E6LL*( Iz/2*(1-cos(ox2)) + Iy/2*(1-cos(ox2)) ); %y nonlinear stiffness including torsion node 2
   
   nonlinFz1=(z1-z2)*E12LLL*( Iy/2*(1-cos(ox1)) + Iz/2*(1-cos(ox1)) ) + -(oy1+oy2)*E6LL*( Iy/2*(1-cos(ox1)) + Iz/2*(1-cos(ox1)) ); %z nonlinear stiffness including torsion node 1
   nonlinFz2=(z2-z1)*E12LLL*( Iy/2*(1-cos(ox2)) + Iz/2*(1-cos(ox2)) ) - -(oy1+oy2)*E6LL*( Iy/2*(1-cos(ox2)) + Iz/2*(1-cos(ox2)) ); %z nonlinear stiffness including torsion node 2
   
   
%   nly=y2*E12LLL*( Iz/2*(1-cos(ox2)) + Iy/2*(1-cos(ox2)) ) - (oz1+oz2)*E6LL*( Iz/2*(1-cos(ox2)) + Iy/2*(1-cos(ox2)) ); 
   
   F=R12'... %rotate force to global
    *[0 
      nonlinFy1
      nonlinFz1
      0
      0
      0
      0
      nonlinFy2
      nonlinFz2
      0
      0
      0];   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   output=[F];
   
   
case 'display'

  displaybeam(q, nodes(1).pos, R, param.l, param.w, param.h);

case 'writebin'

  fid = varargin{1};
  var_ids = varargin{2};
  if (fid < 0)
    output = 1;
  else
    writebeam(fid, nodes(1).pos, R, param.l, param.w, param.h, var_ids);
  end

otherwise

  output = [];

end

