
function [output] = MF_hinge2(flag, R, param, q, t, nodes, varargin);

switch(flag)

case 'vars'

  output.dynamic = {1 {'x' 'y' 'z'};
                    2 {'x' 'y' 'z'}};

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
   
%stiffness matrix, linear model
%   K=a*eye(3);
   K=a*[1 0 0; 0 0 0; 0 0 0];
   R=[R zeros(3); zeros(3) R];
   output = R*[K -K;-K K]*R';


case 'pos'

  % Compute relative positions of beam nodes

  output = R * [0 param.l;
                0 0;
                0 0];            
   
case 'display'
   L=param.l;
   P=R*[L;0;0];
   beta=P(2)+q(5)-q(2);
   gamma=P(3)+q(6)-q(3);
   Q=zeros(12,1);
   Q(1:3,1)=q(1:3);
   Q(7:9,1)=q(4:6);
   Q(4)=0;
   Q(5)=asin(gamma/L)-asin(P(3)/L);
   Q(6)=asin(beta/L)-asin(P(2)/L);
   Q(10)=Q(4);
   Q(11)=Q(5);
   Q(12)=Q(6);
   displaybeam(Q, nodes(1).pos, R, param.l, param.w, param.h);

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

