%Nonlinear beam, stretch
function [output] = MF_beam3dnl(flag, R, param, q, t, nodes, varargin);

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
     J = W*H^3*(16/3-3.36*H/W*(1-(H/W)^4/12));
   else
     J = H*W^3*(16/3-3.36*W/H*(1-(W/H)^4/12));
   end
   % J = (W * H^3 + H*W^3) / 12;       

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
       
   %rotations and displacements
   R3g = R; %Rotate to global
   R3l = R'; %Rotate to local
   dR3l = rot2local(q(4),q(5),q(6)); %dR is the additional rotation of node1 to local. 
   dR3g = dR3l'; %dR is the additional rotation of node1 to global. 
   o = zeros(3,3); %Matrix place-holders. 
   O = zeros(6,6); %Matrix place-holders. 
   Rg = [[R3g o; o R3g] O; O [R3g o; o R3g]]; %12x12 rotation to global. 
   dRg = [[dR3g o; o dR3g] O; O [dR3g o; o dR3g]]; %12x12 rotation to global. 
   Rl = Rg'; %12x12 rotation to local. 
   dRl = dRg'; %12x12 rotation to local. 
   %q is size 24x1, [displacement;velocity] => [q ; qdot]
   qg = q(1:12,1); %global displacement part of q
   qg = [ zeros(3,1); qg(4:6,1); qg(7:9,1)-qg(1:3,1); qg(10:12,1) ]; %global translation to origin 
   ql = Rl * (dRl * qg); %local^2 displacement part of q
%   ql = [ zeros(3,1); ql(4:6,1); ql(7:9,1)-ql(1:3,1); ql(10:12,1) ]; %local translation to origin 

   %geometrical and material parameters
   A = param.w * param.h; %Cross-sectional area of the beam.
   L = param.l; %Beam length, along x-axis. 
   E = param.Youngsmodulus; %Young's modulus.   
       
   %matrix elements   
   Cnl = 0.6*E*A/L^3; %stiffness in x due to y or z displacement.
   Anl = 0.69*E*A/L^3; %stiffness in y due to y or z displacement.
   Bnl = 0.32*E*A/L^2; %moment in z or y due to y or z displacement.
   %Nonlinear stiffness matrix, linear model
   
Bnl=Bnl;   
   
   Kl = [ ... 
   0     Cnl   Cnl   0     0     0,    0    -Cnl  -Cnl   0     0     0;
   0     Anl   Anl   0     0     0,    0    -Anl  -Anl   0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0    -Bnl   0     0     0,    0     0     Bnl   0     0     0;  
   0     Bnl   0     0     0     0,    0    -Bnl   0     0     0     0;   
   0    -Cnl  -Cnl   0     0     0,    0     Cnl   Cnl   0     0     0;
   0    -Anl  -Anl   0     0     0,    0     Anl   Anl   0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0    -Bnl   0     0     0,    0     0     Bnl   0     0     0;
   0     Bnl   0     0     0     0,    0    -Bnl   0     0     0     0];



%   ql = [0 ql(2) ql(3) 0 0 0 0 ql(8) ql(9) 0 0 0]';
   ql = [0 0 0 0 0 0 0 ql(8)-ql(2) ql(9)-ql(3) 0 0 0]';
   Fl = (-1) * Kl * ql.^3;
   Fg = dRg * (Rg * Fl); %global force.
   
   output = Fg;
   
   
%Klocal = Knl;%[K11  K12; K12' K22]; %stiffness in local coordinate frame. 
%   qTranslate2Origin = qGlobal;%[zeros(6,1); qGlobal(7:12,1)-qGlobal(1:6,1)]; %Translate so that node1 is at origin.
%   qTranslate2Origin = qGlobal;%[zeros(3,1); qGlobal(4:6,1); qGlobal(7:9,1)-qGlobal(1:3,1);qGlobal(10:12,1)]; %Translate so that node1 is at origin.
%   qLocalOrigin = ROT * (dROT * qTranslate2Origin); %Rotate so that beam begins along local x.
%   FLocalOrigin = -Klocal * qLocalOrigin.^3; %Forces at local origin.
%   Fglobal = dROT' * (ROT' * FLocalOrigin); %Rotate forces to global frame.
%   output = Fglobal; %return force.   

case 'dFdx'
       
   %rotations and displacements
   R3g = R; %Rotate to global
   R3l = R'; %Rotate to local
   dR3l = rot2local(q(4),q(5),q(6)); %dR is the additional rotation of node1 to local. 
   dR3g = dR3l'; %dR is the additional rotation of node1 to global. 
   o = zeros(3,3); %Matrix place-holders. 
   O = zeros(6,6); %Matrix place-holders. 
   Rg = [[R3g o; o R3g] O; O [R3g o; o R3g]]; %12x12 rotation to global. 
   dRg = [[dR3g o; o dR3g] O; O [dR3g o; o dR3g]]; %12x12 additional rotation to global. 
   Rl = Rg'; %12x12 rotation to local. 
   dRl = dRg'; %12x12 additional rotation to local. 
   %q is size 24x1, [displacement;velocity] => [q ; qdot]
   qg = q(1:12,1); %global displacement part of q
   qg = [ zeros(3,1); qg(4:6,1); qg(7:9,1)-qg(1:3,1); qg(10:12,1) ];%global translation to origin
   ql = Rl * (dRl * qg); %local^2 displacement part of q
%   ql = [ zeros(3,1); ql(4:6,1); ql(7:9,1)-ql(1:3,1); ql(10:12,1) ]; %local translation to origin 
   
   %geometrical and material parameters
   A = param.w * param.h; %Cross-sectional area of the beam.
   L = param.l; %Beam length, along x-axis. 
   E = param.Youngsmodulus; %Young's modulus.   
       
   %matrix elements   
   Cnl = 0.6*E*A/L^3; %stiffness in x due to y or z displacement.
   Anl = 0.69*E*A/L^3; %stiffness in y due to y or z displacement.
   Bnl = 0.32*E*A/L^2; %moment in z or y due to y or z displacement.
   %Nonlinear stiffness matrix, linear model
   
Bnl=Bnl;   
   
   Kl = [ ... 
   0     Cnl   Cnl   0     0     0,    0    -Cnl  -Cnl   0     0     0;
   0     Anl   Anl   0     0     0,    0    -Anl  -Anl   0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0    -Bnl   0     0     0,    0     0     Bnl   0     0     0;  
   0     Bnl   0     0     0     0,    0    -Bnl   0     0     0     0;   
   0    -Cnl  -Cnl   0     0     0,    0     Cnl   Cnl   0     0     0;
   0    -Anl  -Anl   0     0     0,    0     Anl   Anl   0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0    -Bnl   0     0     0,    0     0     Bnl   0     0     0;
   0     Bnl   0     0     0     0,    0    -Bnl   0     0     0     0];


   Fl = (-1) * Kl * ql.^3;
   
%   ql = [0 ql(2) ql(3) 0 0 0 0 ql(8) ql(9) 0 0 0]';
   ql = [0 0 0 0 0 0 0 ql(8)-ql(2) ql(9)-ql(3) 0 0 0]';
   dFdxl = (-1) * 2 * Kl * diag(ql.^2); %local jacobian
   dFdxg = dRg * (Rg * dFdxl * Rl) * dRl; %global jacobian 
   
   output = dFdxg;
  
case 'display'

  displaybeam(q, nodes(1).pos, R, param);
  
otherwise

  output = [];

end

