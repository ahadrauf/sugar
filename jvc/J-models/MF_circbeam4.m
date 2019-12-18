% 3D small-deflection linear beam model
%
% Input parameters:
%   l, w  - beam length and width
%   h - beam height (often specified in process parameters)
%   density, fluid, viscosity, Youngmodulus -
%       material parameters specified as part of the process info
%
% Nodes/variables:
%   The two end nodes are conceptually in the middle of the two
%   w-by-h end faces, which are l units apart.  Each node has
%   the usual spatial displacement variables (x, y, z, rx, ry, rz)

function [output] = MF_circbeam4(flag, R, param, q, t, nodes, varargin);

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
   rotation=R; %store R
   E = param.Youngsmodulus; %Modulus of elasticity.
   H = param.h; %Beam layer thickness, along z-axis.
   W = param.w; %Beam width, along y-axis.
   A = H*W; %Beam crossectional area, yz-plane.
   Nu = param.Poisson; %Poisson's ratio.
   G = E/(2*(1+Nu));	%Shear modulus. GJ = Torsional stiffness 
   if (W>H), J=W*H^3*(16/3-3.36*H/W*(1-(H/W)^4/12))/16; else, J = H*W^3*(16/3-3.36*W/H*(1-(W/H)^4/12))/16; end %The newer version used this different formula for J.  If the two are supposed to be the same, there is a very high relative error (at least ~0.9 for the mirror mode demo). I don't know which one is correct. J2 = 2/7*(W*H)^3/(W^2+H^2);	%Polar second moment of area. abs((J2-J)/J)
   Iy=W*H^3/12;		%Moment about y.
   Iz=H*W^3/12;		%Moment about z.    
   k=zeros(12);
   alpha = 6/5; %form factor for shear. alpha = 6/5 & 10/9 for rectangular and circular cross sections.
   
   AE=A*E; %extensional rigidity
   AG=A*G/alpha; %shear rigidity / alpha. 
   GJ=J*G; %torsional rigidity
   lambda=param.alpha; %curveature range   
   R=param.radius; %radius
   EIz=E*Iz; %flexural rigidity in the xy-plane
   EIy=E*Iy; %flexural rigidity in the xz-plane
   
%stiffness from the xy-plane derivation done in Maple

k(1,1)=...
   2*(-lambda^2*AG*EIz+lambda*alpha*cos(lambda)*sin(lambda)*AE*EIz-lambda*cos(lambda)*sin(lambda)*AG*EIz-alpha*lambda^2*AE*EIz-lambda*R^2*cos(lambda)*sin(lambda)*AE*AG-R^2*lambda^2*AE*AG+2*R^2*AE*AG*sin(lambda)^2)*AE*AG*EIz/(R*(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG));
k(1,2)=...
   -2*(-alpha*lambda*AE*EIz-lambda*AG*EIz+lambda*alpha*cos(lambda)^2*AE*EIz+lambda*cos(lambda)^2*AG*EIz+R^2*lambda*AE*AG-lambda*R^2*cos(lambda)^2*AE*AG-2*R^2*sin(lambda)*AE*AG+2*R^2*cos(lambda)*sin(lambda)*AE*AG)*AE*AG*EIz/(R*(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG));
k(1,6)=...
   -2*(-AG*EIz*sin(lambda)-lambda*alpha*cos(lambda)^2*AE*EIz+lambda*R^2*cos(lambda)*AE*AG+lambda*R^2*cos(lambda)^2*AE*AG-lambda*cos(lambda)^2*AG*EIz-cos(lambda)*sin(lambda)*AG*EIz-3*R^2*cos(lambda)*sin(lambda)*AE*AG+3*R^2*sin(lambda)*AE*AG-2*R^2*lambda*AE*AG+cos(lambda)*alpha*lambda*AE*EIz+cos(lambda)*lambda*AG*EIz-alpha*AE*EIz*sin(lambda)+2*cos(lambda)^2*AG*EIz*sin(lambda)+alpha*cos(lambda)*sin(lambda)*AE*EIz)*AE*AG*EIz/(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG);

k(2,1)=...
   -2*(-alpha*lambda*AE*EIz-lambda*AG*EIz+lambda*alpha*cos(lambda)^2*AE*EIz+lambda*cos(lambda)^2*AG*EIz+R^2*lambda*AE*AG-lambda*R^2*cos(lambda)^2*AE*AG-2*R^2*sin(lambda)*AE*AG+2*R^2*cos(lambda)*sin(lambda)*AE*AG)*AE*AG*EIz/(R*(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-            2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG));
k(2,2)=...      
   2*(-lambda^2*AG*EIz-lambda*alpha*cos(lambda)*sin(lambda)*AE*EIz+lambda*cos(lambda)*sin(lambda)*AG*EIz-alpha*lambda^2*AE*EIz+lambda*R^2*cos(lambda)*sin(lambda)*AE*AG-R^2*lambda^2*AE*AG+2*R^2*AE*AG-4*R^2*cos(lambda)*AE*AG+2*R^2*cos(lambda)^2*AE*AG)*AE*AG*EIz/(R*(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG));
k(2,6)=...   
   -2*(lambda*AG*EIz*sin(lambda)-lambda^2*AG*EIz+alpha*cos(lambda)*sin(lambda)^2*AE*EIz-lambda*alpha*cos(lambda)*sin(lambda)*AE*EIz-cos(lambda)*sin(lambda)^2*AG*EIz+lambda*cos(lambda)*sin(lambda)*AG*EIz+alpha*lambda*AE*EIz*sin(lambda)-alpha*lambda^2*AE*EIz-R^2*cos(lambda)*sin(lambda)^2*AE*AG+lambda*R^2*cos(lambda)*sin(lambda)*AE*AG+lambda*R^2*sin(lambda)*AE*AG-R^2*lambda^2*AE*AG+alpha*AE*EIz-alpha*AE*EIz*cos(lambda)+AG*EIz-AG*EIz*cos(lambda)-alpha*cos(lambda)^2*AE*EIz+alpha*cos(lambda)^3*AE*EIz-cos(lambda)^2*AG*EIz+cos(lambda)^3*AG*EIz+R^2*AE*AG-3*R^2*cos(lambda)*AE*AG+3*R^2*cos(lambda)^2*AE*AG-R^2*cos(lambda)^3*AE*AG)*AE*AG*EIz/(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG);

k(6,1)=...
   -2*(-AG*EIz*sin(lambda)-lambda*alpha*cos(lambda)^2*AE*EIz+lambda*R^2*cos(lambda)*AE*AG+lambda*R^2*cos(lambda)^2*AE*AG-lambda*cos(lambda)^2*AG*EIz-cos(lambda)*sin(lambda)*AG*EIz-3*R^2*cos(lambda)*sin(lambda)*AE*AG+3*R^2*sin(lambda)*AE*AG-2*R^2*lambda*AE*AG+cos(lambda)*alpha*lambda*AE*EIz+cos(lambda)*lambda*AG*EIz-alpha*AE*EIz*sin(lambda)+2*cos(lambda)^2*AG*EIz*sin(lambda)+alpha*cos(lambda)*sin(lambda)*AE*EIz)*AE*AG*EIz/(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG);
k(6,2)=...   
   -2*(lambda*AG*EIz*sin(lambda)-lambda^2*AG*EIz+alpha*cos(lambda)*sin(lambda)^2*AE*EIz-lambda*alpha*cos(lambda)*sin(lambda)*AE*EIz-cos(lambda)*sin(lambda)^2*AG*EIz+lambda*cos(lambda)*sin(lambda)*AG*EIz+alpha*lambda*AE*EIz*sin(lambda)-alpha*lambda^2*AE*EIz-R^2*cos(lambda)*sin(lambda)^2*AE*AG+lambda*R^2*cos(lambda)*sin(lambda)*AE*AG+lambda*R^2*sin(lambda)*AE*AG-R^2*lambda^2*AE*AG+alpha*AE*EIz-alpha*AE*EIz*cos(lambda)+AG*EIz-AG*EIz*cos(lambda)-alpha*cos(lambda)^2*AE*EIz+alpha*cos(lambda)^3*AE*EIz-cos(lambda)^2*AG*EIz+cos(lambda)^3*AG*EIz+R^2*AE*AG-3*R^2*cos(lambda)*AE*AG+3*R^2*cos(lambda)^2*AE*AG-R^2*cos(lambda)^3*AE*AG)*AE*AG*EIz/(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG);
k(6,6)=...   
   (4*alpha*lambda*AE^2*EIz*R^2*sin(lambda)*AG+4*alpha*cos(lambda)^3*AE^2*EIz*R^2*AG+2*alpha*AE^2*EIz*R^2*AG+AG^2*EIz^2-4*alpha*lambda^2*AE^2*EIz*R^2*AG-3*R^4*lambda^2*AE^2*AG^2-4*lambda^2*AG^2*EIz*R^2*AE-2*lambda^2*AG*EIz^2*alpha*AE-alpha^2*lambda^2*AE^2*EIz^2+2*alpha*AE*EIz^2*AG+2*AG^2*EIz*R^2*AE+alpha^2*AE^2*EIz^2-lambda^2*AG^2*EIz^2+R^4*AE^2*AG^2-2*AG^2*EIz^2*cos(lambda)^2+cos(lambda)^4*AG^2*EIz^2+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*R^2*AE-4*R^4*cos(lambda)*sin(lambda)^2*AE^2*AG^2-4*alpha*AE^2*EIz*R^2*cos(lambda)*AG-4*alpha*AE*EIz^2*cos(lambda)^2*AG+2*alpha*cos(lambda)^4*AE*EIz^2*AG-2*alpha*cos(lambda)^4*AE^2*EIz*R^2*AG+4*cos(lambda)^3*AG^2*EIz*R^2*AE+4*lambda*AG^2*EIz*R^2*sin(lambda)*AE+2*lambda*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2+2*R^4*cos(lambda)*sin(lambda)*AE^2*AG^2*lambda+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2+4*R^4*lambda*AE^2*AG^2*sin(lambda)-4*AG^2*EIz*R^2*cos(lambda)*AE-2*cos(lambda)^4*AG^2*EIz*R^2*AE-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*AG+4*alpha*cos(lambda)*sin(lambda)^2*AE^2*EIz*R^2*AG-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*R^2*AG-2*alpha*cos(lambda)*sin(lambda)*AE^2*EIz*R^2*lambda*AG-4*cos(lambda)*sin(lambda)^2*AG^2*EIz*R^2*AE-2*alpha^2*AE^2*EIz^2*cos(lambda)^2+alpha^2*cos(lambda)^4*AE^2*EIz^2-4*R^4*AE^2*AG^2*cos(lambda)+6*R^4*AE^2*AG^2*cos(lambda)^2-4*R^4*cos(lambda)^3*AE^2*AG^2+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2+R^4*cos(lambda)^4*AE^2*AG^2)*EIz/(R*(-2*alpha*lambda^3*AE^2*EIz*R^2*AG-2*lambda^3*AG*EIz^2*alpha*AE-2*lambda^3*AG^2*EIz*R^2*AE-lambda^3*AG^2*EIz^2+AG^2*EIz^2*lambda-alpha^2*lambda^3*AE^2*EIz^2-R^4*lambda^3*AE^2*AG^2+alpha^2*AE^2*EIz^2*lambda+3*R^4*AE^2*AG^2*lambda-2*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+2*alpha*cos(lambda)^3*AE^2*EIz*R^2*sin(lambda)*AG-2*cos(lambda)^4*AG^2*EIz*lambda*R^2*AE-8*cos(lambda)^2*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)^4*AE*EIz^2*lambda*AG+6*cos(lambda)^3*AG^2*EIz*R^2*sin(lambda)*AE+R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda-4*alpha*AE*EIz^2*lambda*cos(lambda)^2*AG+cos(lambda)^2*sin(lambda)^2*AG^2*EIz^2*lambda-2*AG^2*EIz^2*lambda*cos(lambda)^2+2*alpha*lambda*AE^2*EIz*R^2*AG*sin(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG+2*lambda*AG^2*EIz*R^2*AE*sin(lambda)^2+2*R^4*lambda*AE^2*AG^2*sin(lambda)^2-6*alpha*AE^2*EIz*R^2*cos(lambda)*sin(lambda)*AG+4*AG^2*EIz*R^2*sin(lambda)*AE+2*alpha*cos(lambda)*sin(lambda)^3*AE^2*EIz*R^2*AG-2*cos(lambda)*sin(lambda)^3*AG^2*EIz*R^2*AE+4*alpha*AE^2*EIz*R^2*sin(lambda)*AG+6*AG^2*EIz*lambda*R^2*cos(lambda)^2*AE+R^4*cos(lambda)^4*AE^2*AG^2*lambda-2*R^4*cos(lambda)^3*AE^2*AG^2*sin(lambda)+6*R^4*AE^2*AG^2*cos(lambda)*sin(lambda)-4*R^2*cos(lambda)*AE^2*AG*alpha*lambda*EIz-4*R^2*cos(lambda)*AE*AG^2*lambda*EIz-4*R^4*cos(lambda)*AE^2*AG^2*lambda-2*R^4*cos(lambda)*sin(lambda)^3*AE^2*AG^2+6*alpha*AE^2*EIz*lambda*R^2*cos(lambda)^2*AG+alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda+cos(lambda)^4*AG^2*EIz^2*lambda+2*cos(lambda)^2*sin(lambda)^2*AG^2*EIz*lambda*R^2*AE-2*AG^2*EIz*R^2*cos(lambda)*sin(lambda)*AE+alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda-4*R^4*AE^2*AG^2*sin(lambda)-2*alpha^2*AE^2*EIz^2*lambda*cos(lambda)^2-2*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+2*alpha*AE*EIz^2*lambda*AG));

%stiffness from the xz-plane derivation done in Maple
k(3,3)=...
   (-2*lambda^2*GJ*EIy-lambda^2*GJ^2-lambda^2*EIy^2-2*cos(lambda)^2*sin(lambda)^2*GJ*EIy+cos(lambda)^2*sin(lambda)^2*GJ^2+cos(lambda)^2*sin(lambda)^2*EIy^2+GJ^2+2*EIy*GJ-2*cos(lambda)^2*GJ^2-4*GJ*cos(lambda)^2*EIy+EIy^2-2*cos(lambda)^2*EIy^2+cos(lambda)^4*GJ^2+2*cos(lambda)^4*GJ*EIy+cos(lambda)^4*EIy^2)*GJ*AG/(R*(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2));
k(3,4)=...
   -2*EIy*(2*cos(lambda)^3*sin(lambda)*GJ-2*cos(lambda)^2*sin(lambda)*GJ+cos(lambda)*lambda*GJ+sin(lambda)*GJ-cos(lambda)*sin(lambda)*GJ-lambda*GJ+cos(lambda)*EIy*lambda-cos(lambda)*sin(lambda)*EIy-lambda*EIy+sin(lambda)*EIy)*GJ*AG/(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-      10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2);
k(3,5)=...
   -(2*EIy*GJ+EIy^2-2*GJ*cos(lambda)^2*EIy+GJ^2-2*cos(lambda)^2*GJ^2-2*cos(lambda)*EIy*GJ-2*cos(lambda)*EIy^2+2*cos(lambda)^3*EIy*GJ+2*cos(lambda)^3*EIy^2-cos(lambda)^4*EIy^2+cos(lambda)^4*GJ^2+lambda^2*EIy^2+4*lambda*GJ*cos(lambda)*sin(lambda)*EIy-lambda^2*GJ^2-cos(lambda)^2*sin(lambda)^2*EIy^2-2*sin(lambda)*EIy*lambda*GJ-2*sin(lambda)*EIy^2*lambda-2*sin(lambda)^2*EIy*cos(lambda)*GJ+2*sin(lambda)^2*EIy^2*cos(lambda)+cos(lambda)^2*sin(lambda)^2*GJ^2)*GJ*AG/(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2);

k(4,3)=...
   -2*EIy*(2*cos(lambda)^3*sin(lambda)*GJ-2*cos(lambda)^2*sin(lambda)*GJ+cos(lambda)*lambda*GJ+sin(lambda)*GJ-cos(lambda)*sin(lambda)*GJ-lambda*GJ+cos(lambda)*EIy*lambda-cos(lambda)*sin(lambda)*EIy-lambda*EIy+sin(lambda)*EIy)*GJ*AG/(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2);
k(4,4)=...
   2*EIy*(alpha*lambda*GJ^2*cos(lambda)*sin(lambda)-alpha*lambda*GJ*EIy*cos(lambda)*sin(lambda)-alpha*lambda^2*GJ*EIy+R^2*lambda*AG*cos(lambda)*sin(lambda)*GJ-R^2*lambda^2*EIy*AG-3*R^2*lambda^2*AG*GJ+2*R^2*AG*sin(lambda)^2*EIy+2*R^2*sin(lambda)^2*cos(lambda)^2*GJ*AG-4*R^2*sin(lambda)^2*AG*cos(lambda)*GJ+4*R^2*sin(lambda)*AG*lambda*GJ-R^2*lambda*EIy*AG*cos(lambda)*sin(lambda)-alpha*lambda^2*GJ^2)*GJ/(R*(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2));
k(4,5)=...
   2*EIy*(-R^2*lambda*EIy*AG*cos(lambda)^2+2*R^2*lambda*EIy*AG-R^2*EIy*AG*cos(lambda)*lambda+R^2*AG*cos(lambda)*lambda*GJ+3*R^2*sin(lambda)*cos(lambda)*EIy*AG+R^2*sin(lambda)*cos(lambda)*GJ*AG-alpha*lambda*GJ^2*cos(lambda)^2-alpha*lambda*GJ*EIy*cos(lambda)^2+alpha*lambda*GJ^2+alpha*lambda*GJ*EIy-3*R^2*sin(lambda)*GJ*AG+2*R^2*sin(lambda)*AG*GJ*cos(lambda)^2-3*R^2*sin(lambda)*EIy*AG-3*R^2*lambda*AG*GJ*cos(lambda)^2+2*R^2*lambda*GJ*AG)*GJ/(R*(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2));

k(5,3)=...
   -(2*EIy*GJ+EIy^2-2*GJ*cos(lambda)^2*EIy+GJ^2-2*cos(lambda)^2*GJ^2-2*cos(lambda)*EIy*GJ-2*cos(lambda)*EIy^2+2*cos(lambda)^3*EIy*GJ+2*cos(lambda)^3*EIy^2-cos(lambda)^4*EIy^2+cos(lambda)^4*GJ^2+lambda^2*EIy^2+4*lambda*GJ*cos(lambda)*sin(lambda)*EIy-lambda^2*GJ^2-cos(lambda)^2*sin(lambda)^2*EIy^2-2*sin(lambda)*EIy*lambda*GJ-2*sin(lambda)*EIy^2*lambda-2*sin(lambda)^2*EIy*cos(lambda)*GJ+2*sin(lambda)^2*EIy^2*cos(lambda)+cos(lambda)^2*sin(lambda)^2*GJ^2)*GJ*AG/(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2);
k(5,4)=...
   2*EIy*(-R^2*lambda*EIy*AG*cos(lambda)^2+2*R^2*lambda*EIy*AG-R^2*EIy*AG*cos(lambda)*lambda+R^2*AG*cos(lambda)*lambda*GJ+3*R^2*sin(lambda)*cos(lambda)*EIy*AG+R^2*sin(lambda)*cos(lambda)*GJ*AG-alpha*lambda*GJ^2*cos(lambda)^2-alpha*lambda*GJ*EIy*cos(lambda)^2+alpha*lambda*GJ^2+alpha*lambda*GJ*EIy-3*R^2*sin(lambda)*GJ*AG+2*R^2*sin(lambda)*AG*GJ*cos(lambda)^2-3*R^2*sin(lambda)*EIy*AG-3*R^2*lambda*AG*GJ*cos(lambda)^2+2*R^2*lambda*GJ*AG)*GJ/(R*(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2));
k(5,5)=...
   (-2*alpha*lambda^2*GJ*EIy^2-R^2*lambda^2*GJ^2*AG-3*R^2*lambda^2*EIy^2*AG-4*R^2*lambda^2*EIy*AG*GJ-2*alpha*lambda^2*GJ^2*EIy+2*R^2*EIy*AG*GJ+R^2*AG*GJ^2-4*R^2*AG*cos(lambda)*EIy^2-4*R^2*AG*cos(lambda)^3*EIy^2+R^2*AG*EIy^2+R^2*AG*cos(lambda)^4*GJ^2+4*R^2*AG*cos(lambda)^3*EIy*GJ-2*R^2*AG*cos(lambda)^4*GJ*EIy-4*R^2*AG*cos(lambda)*EIy*GJ+6*R^2*AG*cos(lambda)^2*EIy^2+2*R^2*lambda*EIy^2*AG*cos(lambda)*sin(lambda)+R^2*AG*cos(lambda)^4*EIy^2-4*R^2*sin(lambda)^2*EIy^2*AG*cos(lambda)+4*R^2*sin(lambda)*EIy*AG*lambda*GJ+4*R^2*sin(lambda)*EIy^2*AG*lambda+R^2*sin(lambda)^2*cos(lambda)^2*GJ^2*AG-2*R^2*lambda*EIy*AG*cos(lambda)*sin(lambda)*GJ-2*R^2*sin(lambda)^2*cos(lambda)^2*GJ*AG*EIy+2*alpha*lambda*GJ*EIy^2*cos(lambda)*sin(lambda)+4*R^2*sin(lambda)^2*EIy*AG*cos(lambda)*GJ-2*alpha*lambda*GJ^2*EIy*cos(lambda)*sin(lambda)-2*R^2*AG*cos(lambda)^2*GJ^2+R^2*sin(lambda)^2*cos(lambda)^2*EIy^2*AG)*GJ/(R*(2*R^2*AG*sin(lambda)^2*EIy^2*lambda+3*R^2*lambda*EIy^2*AG-4*R^2*sin(lambda)*AG*cos(lambda)^4*GJ^2-2*R^2*AG*sin(lambda)^3*EIy^2*cos(lambda)-4*R^2*EIy^2*AG*cos(lambda)*lambda-2*alpha*lambda^3*GJ^2*EIy-alpha*lambda^3*GJ^3-alpha*lambda^3*GJ*EIy^2+alpha*lambda*GJ^3+2*alpha*lambda*GJ^2*EIy+alpha*lambda*GJ*EIy^2-4*R^2*lambda^3*EIy*AG*GJ-3*R^2*lambda^3*AG*GJ^2+3*R^2*lambda*AG*GJ^2-R^2*lambda^3*EIy^2*AG+6*R^2*lambda*EIy*AG*GJ-2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG*lambda^2+4*R^2*lambda*EIy*AG*cos(lambda)^4*GJ+3*R^2*lambda*AG*cos(lambda)^4*GJ^2-6*R^2*lambda*EIy*AG*GJ*cos(lambda)^2-6*R^2*lambda*AG*cos(lambda)^2*GJ^2+R^2*lambda*EIy^2*AG*cos(lambda)^2*sin(lambda)^2-4*R^2*sin(lambda)*EIy^2*AG+8*R^2*sin(lambda)*EIy*AG*GJ*cos(lambda)^2+8*R^2*sin(lambda)*AG*cos(lambda)^2*GJ^2+2*alpha*lambda*GJ^2*EIy*cos(lambda)^4+4*R^2*sin(lambda)^3*EIy*AG*cos(lambda)^2*GJ-8*R^2*sin(lambda)*EIy*AG*GJ-4*R^2*sin(lambda)*AG*GJ^2+4*R^2*sin(lambda)*AG*lambda^2*GJ^2+4*R^2*sin(lambda)*EIy*AG*lambda^2*GJ+alpha*lambda*GJ*EIy^2*cos(lambda)^4+alpha*lambda*GJ^3*cos(lambda)^4-2*alpha*lambda*GJ*EIy^2*cos(lambda)^2+alpha*lambda*GJ*EIy^2*cos(lambda)^2*sin(lambda)^2+alpha*lambda*GJ^3*cos(lambda)^2*sin(lambda)^2-4*alpha*lambda*GJ^2*EIy*cos(lambda)^2+2*R^2*sin(lambda)*cos(lambda)^5*GJ^2*AG-2*alpha*lambda*GJ^2*EIy*cos(lambda)^2*sin(lambda)^2-2*alpha*lambda*GJ^3*cos(lambda)^2+R^2*lambda*EIy^2*AG*cos(lambda)^4-2*R^2*sin(lambda)*cos(lambda)^5*GJ*AG*EIy-10*R^2*sin(lambda)*cos(lambda)^3*GJ*AG*EIy+8*R^2*sin(lambda)*cos(lambda)*GJ*AG*EIy-4*R^2*sin(lambda)*cos(lambda)^3*GJ^2*AG-2*R^2*sin(lambda)*cos(lambda)^3*EIy^2*AG+6*R^2*sin(lambda)*cos(lambda)*EIy^2*AG+2*R^2*sin(lambda)*cos(lambda)*GJ^2*AG-2*R^2*sin(lambda)^3*cos(lambda)^3*GJ*AG*EIy+4*R^2*sin(lambda)*EIy*AG*cos(lambda)^4*GJ+2*R^2*sin(lambda)^3*cos(lambda)^3*GJ^2*AG-4*R^2*EIy*AG*cos(lambda)*lambda*GJ+2*R^2*AG*sin(lambda)^2*EIy*lambda*GJ+2*R^2*AG*sin(lambda)^3*EIy*cos(lambda)*GJ+3*R^2*lambda*AG*cos(lambda)^2*sin(lambda)^2*GJ^2+2*R^2*sin(lambda)*cos(lambda)*GJ*AG*lambda^2*EIy-8*R^2*AG*lambda*EIy*sin(lambda)^2*cos(lambda)*GJ-4*R^2*sin(lambda)^3*AG*cos(lambda)^2*GJ^2));

%Other xy-plane stiffness coefficients from equilibrium laws
%   k([1,2,6],[1,2,6])=Kxy;
   k(7,1)=-k(1,1); 
   k(8,1)=-k(2,1); 
   k(7,2)=-k(2,1);
   k(8,2)=-k(2,2);
   k(7,6)=-k(6,1);
   k(8,6)=-k(6,2);   
   k(7,7)=-k(7,1); 
   k(8,7)=-k(7,2);
   k(8,8)=-k(8,2);
   k(8,7)=-k(7,1); 	
   k(12,8)=-k(12,2); 
   k(12,7)=-k(12,1); 
   k(12,1)=-k(1,1)*R*sin(lambda)-k(1,2)*R*(1-cos(lambda))-k(1,6);
   k(12,2)=-k(2,1)*R*sin(lambda)-k(2,2)*R*(1-cos(lambda))-k(2,6);
   k(12,6)=-k(6,1)*R*sin(lambda)-k(6,2)*R*(1-cos(lambda))-k(6,6);
   k(12,7)=-k(7,1)*R*sin(lambda)-k(7,2)*R*(1-cos(lambda))-k(7,6);
   k(12,8)=-k(8,1)*R*sin(lambda)-k(8,2)*R*(1-cos(lambda))-k(8,6);
   k(12,12)=-k(12,1)*R*sin(lambda)-k(12,2)*R*(1-cos(lambda))-k(12,6); 	
   
   %Other xz-plane stiffness coefficients from equilibrium laws
%   k([3,4,5],[3,4,5])=Kxz;
   k(10,3)=-k(4,3);
   k(11,3)=-k(5,3);
   k(10,4)=-k(4,4);
   k(11,4)=-k(5,4);
   k(10,5)=-k(5,4);
   k(11,5)=-k(5,5);
   k(10,10)=k(4,4);  	
   k(11,10)=k(5,4);
   k(11,11)=k(5,5);  	
   k(9,3)=-k(3,4)*R*sin(lambda)-k(3,5)*R*(1-cos(lambda))-k(3,3); 
   k(9,4)=-k(4,4)*R*sin(lambda)-k(4,5)*R*(1-cos(lambda))-k(4,3); 
   k(9,5)=-k(5,4)*R*sin(lambda)-k(5,5)*R*(1-cos(lambda))-k(5,3); 
   k(9,9)=-k(9,4)*R*sin(lambda)-k(9,5)*R*(1-cos(lambda))-k(9,3);  	
   k(9,10)=-k(10,4)*R*sin(lambda)-k(10,5)*R*(1-cos(lambda))-k(10,3); 
   k(9,11)=-k(11,4)*R*sin(lambda)-k(11,5)*R*(1-cos(lambda))-k(11,3); 
   
   %fill in the symmetric part
   for i=1:12
      for j=1:12
         if k(i,j)~=0,k(j,i)=k(i,j);end
      end
   end
   
   %rotation matrix
R=rotation;
o3=zeros(3);
R=[R,o3;o3,R];
o6=zeros(6);
R=[R,o6;o6,R];
output = R*k*R';

case 'pos'

  % Compute relative positions of beam nodes

  output = R * [0 param.radius*sin(param.alpha);
                0 param.radius*(1-cos(param.alpha));
                0 0];            

case 'F'

   output = [];
   
case 'display'

   displaycircbeam3(q, nodes(1).pos, R, param.w, param.h, param.radius, param.alpha);

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

