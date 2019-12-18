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

function [output] = MF_circbeam6(flag, R, param, q, t, nodes, varargin);

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
k(1,1)=36*(R^2*lambda^2*AE*AG+alpha*lambda^2*AE*EIz+lambda*alpha*cos(lambda)*sin(lambda)*AE*EIz-lambda*R^2*cos(lambda)*sin(lambda)*AE*AG-lambda*cos(lambda)*sin(lambda)*EIz*AG+lambda^2*EIz*AG-2*R^2*AE*AG*cos(lambda)^2+4*R^2*AE*AG*cos(lambda)-2*R^2*AE*AG)*AE*EIz*AG/(R*(-18*EIz^2*AG^2*lambda+36*alpha*lambda^3*AE*EIz^2*AG+36*R^2*lambda^3*AE*AG^2*EIz+36*R^2*lambda^3*AE^2*AG*alpha*EIz+18*R^4*lambda^3*AE^2*AG^2-54*R^4*lambda*AE^2*AG^2+18*alpha^2*lambda^3*AE^2*EIz^2-18*alpha^2*AE^2*EIz^2*lambda+18*lambda^3*EIz^2*AG^2-24*R^2*cos(lambda)^2*sin(lambda)*AE*AG^2*lambda^2*EIz-36*R^4*cos(lambda)^4*sin(lambda)*AE^2*AG^2-12*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2+24*cos(lambda)^5*EIz*AG^2*lambda*R^2*AE-72*cos(lambda)^2*EIz^2*AG*alpha*lambda*AE-36*cos(lambda)^4*EIz*AG^2*lambda*R^2*AE-48*cos(lambda)^3*EIz*AG^2*lambda*R^2*AE-45*cos(lambda)^2*sin(lambda)^2*EIz*AG^2*lambda*R^2*AE+36*cos(lambda)^3*sin(lambda)*EIz*AG^2*R^2*AE-108*cos(lambda)*sin(lambda)*EIz*AG^2*R^2*AE-11*R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda+36*lambda*EIz*AG^2*R^2*AE*cos(lambda)^2+96*lambda*EIz*AG^2*R^2*AE*cos(lambda)-9*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2*lambda^2-9*R^2*cos(lambda)^3*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-9*R^2*cos(lambda)^4*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz-7*R^4*cos(lambda)^4*sin(lambda)^2*AE^2*AG^2*lambda+9*R^2*cos(lambda)^4*sin(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*cos(lambda)^3*sin(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*cos(lambda)^5*sin(lambda)*AE^2*AG^2-9*R^4*sin(lambda)^3*cos(lambda)*AE^2*AG^2*lambda^2-9*R^2*sin(lambda)^3*cos(lambda)*AE^2*AG*alpha*lambda^2*EIz-24*R^2*sin(lambda)^4*AE^2*AG*lambda*alpha*cos(lambda)*EIz+24*R^4*sin(lambda)^4*AE^2*AG^2*lambda*cos(lambda)+24*R^2*sin(lambda)^4*AE*AG^2*lambda*cos(lambda)*EIz+48*R^2*cos(lambda)^3*sin(lambda)^2*AE*AG^2*lambda*EIz+45*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+36*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG-36*alpha*cos(lambda)^3*sin(lambda)*AE^2*EIz*R^2*AG+108*alpha*cos(lambda)*sin(lambda)*AE^2*EIz*R^2*AG-24*R^4*sin(lambda)^3*AE^2*AG^2*lambda^2+96*R^4*lambda*AE^2*AG^2*cos(lambda)-108*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)^2+48*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)-24*R^2*sin(lambda)^3*AE*AG^2*lambda^2*EIz-36*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)^2+6*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)+24*R^4*AE^2*AG^2*sin(lambda)*lambda^2+24*R^2*AE^2*AG*sin(lambda)*alpha*lambda^2*EIz+48*R^2*AE^2*AG*sin(lambda)^2*lambda*alpha*cos(lambda)*EIz-48*R^4*AE^2*AG^2*sin(lambda)^2*lambda*cos(lambda)-48*R^2*AE*AG^2*sin(lambda)^2*lambda*cos(lambda)*EIz+24*R^2*AE*AG^2*sin(lambda)*lambda^2*EIz+84*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)^2-78*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)-24*R^4*cos(lambda)^2*sin(lambda)*AE^2*AG^2*lambda^2-24*R^2*cos(lambda)^2*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-48*R^2*cos(lambda)^3*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz+48*R^4*cos(lambda)^3*sin(lambda)^2*AE^2*AG^2*lambda-2*R^4*AE^2*AG^2*cos(lambda)^4*lambda-48*R^4*AE^2*AG^2*cos(lambda)^3*lambda+9*R^2*lambda^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)*EIz+9*R^4*lambda^2*AE^2*AG^2*cos(lambda)*sin(lambda)+9*R^2*lambda^2*AE*AG^2*cos(lambda)*sin(lambda)*EIz+72*EIz*AG^2*R^2*AE*sin(lambda)-72*alpha*AE^2*EIz*R^2*AG*sin(lambda)-18*alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda+36*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+36*alpha^2*cos(lambda)^2*AE^2*EIz^2*lambda+48*R^2*AE^2*AG*cos(lambda)^3*lambda*alpha*EIz-8*R^4*cos(lambda)^6*AE^2*AG^2*lambda-24*R^2*cos(lambda)^5*AE^2*AG*lambda*alpha*EIz+24*R^4*cos(lambda)^5*AE^2*AG^2*lambda+R^4*sin(lambda)^4*cos(lambda)^2*AE^2*AG^2*lambda+9*R^2*sin(lambda)^4*cos(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*sin(lambda)^3*cos(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*sin(lambda)^3*cos(lambda)^3*AE^2*AG^2-18*cos(lambda)^2*sin(lambda)^2*EIz^2*AG^2*lambda-8*R^4*lambda*AE^2*AG^2*cos(lambda)^2-18*alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda-24*R^2*sin(lambda)^3*AE^2*AG*alpha*lambda^2*EIz-9*R^2*sin(lambda)^4*cos(lambda)^2*AE^2*AG*lambda*alpha*EIz-36*R^2*AE*AG^2*sin(lambda)^2*lambda*EIz-36*R^4*AE^2*AG^2*sin(lambda)^2*lambda-...
   36*R^2*AE^2*AG*sin(lambda)^2*alpha*lambda*EIz-36*R^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)^3*EIz+36*R^2*AE*AG^2*cos(lambda)*sin(lambda)^3*EIz+48*R^4*sin(lambda)^3*AE^2*AG^2+24*R^4*AE^2*AG^2*sin(lambda)-18*cos(lambda)^4*EIz^2*AG^2*lambda+36*cos(lambda)^2*EIz^2*AG^2*lambda+36*EIz^2*AG*alpha*lambda*AE-72*lambda*EIz*AG^2*R^2*AE+36*cos(lambda)^4*EIz^2*AG*lambda*alpha*AE));
k(1,2)=-12*(-3*lambda*cos(lambda)^2*EIz*AG-2*lambda*R^2*AE*AG*cos(lambda)+2*lambda*R^2*cos(lambda)^3*AE*AG+3*lambda*alpha*cos(lambda)^2*AE*EIz-3*lambda*R^2*AE*AG*cos(lambda)^2+2*lambda*R^2*sin(lambda)^2*cos(lambda)*AE*AG+3*R^2*lambda*AE*AG+3*lambda*EIz*AG-3*alpha*lambda*AE*EIz+6*R^2*cos(lambda)*sin(lambda)*AE*AG-6*R^2*AE*AG*sin(lambda))*AE*EIz*AG/(R*(-18*EIz^2*AG^2*lambda+36*alpha*lambda^3*AE*EIz^2*AG+36*R^2*lambda^3*AE*AG^2*EIz+36*R^2*lambda^3*AE^2*AG*alpha*EIz+18*R^4*lambda^3*AE^2*AG^2-54*R^4*lambda*AE^2*AG^2+18*alpha^2*lambda^3*AE^2*EIz^2-18*alpha^2*AE^2*EIz^2*lambda+18*lambda^3*EIz^2*AG^2-24*R^2*cos(lambda)^2*sin(lambda)*AE*AG^2*lambda^2*EIz-36*R^4*cos(lambda)^4*sin(lambda)*AE^2*AG^2-12*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2+24*cos(lambda)^5*EIz*AG^2*lambda*R^2*AE-72*cos(lambda)^2*EIz^2*AG*alpha*lambda*AE-36*cos(lambda)^4*EIz*AG^2*lambda*R^2*AE-48*cos(lambda)^3*EIz*AG^2*lambda*R^2*AE-45*cos(lambda)^2*sin(lambda)^2*EIz*AG^2*lambda*R^2*AE+36*cos(lambda)^3*sin(lambda)*EIz*AG^2*R^2*AE-108*cos(lambda)*sin(lambda)*EIz*AG^2*R^2*AE-11*R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda+36*lambda*EIz*AG^2*R^2*AE*cos(lambda)^2+96*lambda*EIz*AG^2*R^2*AE*cos(lambda)-9*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2*lambda^2-9*R^2*cos(lambda)^3*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-9*R^2*cos(lambda)^4*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz-7*R^4*cos(lambda)^4*sin(lambda)^2*AE^2*AG^2*lambda+9*R^2*cos(lambda)^4*sin(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*cos(lambda)^3*sin(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*cos(lambda)^5*sin(lambda)*AE^2*AG^2-9*R^4*sin(lambda)^3*cos(lambda)*AE^2*AG^2*lambda^2-9*R^2*sin(lambda)^3*cos(lambda)*AE^2*AG*alpha*lambda^2*EIz-24*R^2*sin(lambda)^4*AE^2*AG*lambda*alpha*cos(lambda)*EIz+24*R^4*sin(lambda)^4*AE^2*AG^2*lambda*cos(lambda)+24*R^2*sin(lambda)^4*AE*AG^2*lambda*cos(lambda)*EIz+48*R^2*cos(lambda)^3*sin(lambda)^2*AE*AG^2*lambda*EIz+45*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+36*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG-36*alpha*cos(lambda)^3*sin(lambda)*AE^2*EIz*R^2*AG+108*alpha*cos(lambda)*sin(lambda)*AE^2*EIz*R^2*AG-24*R^4*sin(lambda)^3*AE^2*AG^2*lambda^2+96*R^4*lambda*AE^2*AG^2*cos(lambda)-108*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)^2+48*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)-24*R^2*sin(lambda)^3*AE*AG^2*lambda^2*EIz-36*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)^2+6*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)+24*R^4*AE^2*AG^2*sin(lambda)*lambda^2+24*R^2*AE^2*AG*sin(lambda)*alpha*lambda^2*EIz+48*R^2*AE^2*AG*sin(lambda)^2*lambda*alpha*cos(lambda)*EIz-48*R^4*AE^2*AG^2*sin(lambda)^2*lambda*cos(lambda)-48*R^2*AE*AG^2*sin(lambda)^2*lambda*cos(lambda)*EIz+24*R^2*AE*AG^2*sin(lambda)*lambda^2*EIz+84*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)^2-78*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)-24*R^4*cos(lambda)^2*sin(lambda)*AE^2*AG^2*lambda^2-24*R^2*cos(lambda)^2*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-48*R^2*cos(lambda)^3*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz+48*R^4*cos(lambda)^3*sin(lambda)^2*AE^2*AG^2*lambda-2*R^4*AE^2*AG^2*cos(lambda)^4*lambda-48*R^4*AE^2*AG^2*cos(lambda)^3*lambda+9*R^2*lambda^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)*EIz+9*R^4*lambda^2*AE^2*AG^2*cos(lambda)*sin(lambda)+9*R^2*lambda^2*AE*AG^2*cos(lambda)*sin(lambda)*EIz+72*EIz*AG^2*R^2*AE*sin(lambda)-72*alpha*AE^2*EIz*R^2*AG*sin(lambda)-18*alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda+36*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+36*alpha^2*cos(lambda)^2*AE^2*EIz^2*lambda+48*R^2*AE^2*AG*cos(lambda)^3*lambda*alpha*EIz-8*R^4*cos(lambda)^6*AE^2*AG^2*lambda-24*R^2*cos(lambda)^5*AE^2*AG*lambda*alpha*EIz+24*R^4*cos(lambda)^5*AE^2*AG^2*lambda+R^4*sin(lambda)^4*cos(lambda)^2*AE^2*AG^2*lambda+9*R^2*sin(lambda)^4*cos(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*sin(lambda)^3*cos(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*sin(lambda)^3*cos(lambda)^3*AE^2*AG^2-18*cos(lambda)^2*sin(lambda)^2*EIz^2*AG^2*lambda-8*R^4*lambda*AE^2*AG^2*cos(lambda)^2-18*alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda-24*R^2*sin(lambda)^3*AE^2*AG*alpha*lambda^2*EIz-9*R^2*sin(lambda)^4*cos(lambda)^2*AE^2*AG*lambda*alpha*EIz-...
   36*R^2*AE*AG^2*sin(lambda)^2*lambda*EIz-36*R^4*AE^2*AG^2*sin(lambda)^2*lambda-36*R^2*AE^2*AG*sin(lambda)^2*alpha*lambda*EIz-36*R^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)^3*EIz+36*R^2*AE*AG^2*cos(lambda)*sin(lambda)^3*EIz+48*R^4*sin(lambda)^3*AE^2*AG^2+24*R^4*AE^2*AG^2*sin(lambda)-18*cos(lambda)^4*EIz^2*AG^2*lambda+36*cos(lambda)^2*EIz^2*AG^2*lambda+36*EIz^2*AG*alpha*lambda*AE-72*lambda*EIz*AG^2*R^2*AE+36*cos(lambda)^4*EIz^2*AG*lambda*alpha*AE));
k(1,6)=-12*(3*EIz*AG-3*alpha*AE*EIz+3*R^2*lambda^2*AE*AG+3*alpha*lambda^2*AE*EIz-3*R^2*AE*AG-7*R^2*AE*AG*cos(lambda)^2+7*R^2*AE*AG*cos(lambda)+3*lambda*alpha*cos(lambda)*sin(lambda)*AE*EIz-3*lambda*R^2*cos(lambda)*sin(lambda)*AE*AG-3*lambda*cos(lambda)*sin(lambda)*EIz*AG+3*lambda^2*EIz*AG-2*R^2*sin(lambda)^2*cos(lambda)^2*AE*AG-3*alpha*cos(lambda)*sin(lambda)^2*AE*EIz-3*sin(lambda)*alpha*lambda*AE*EIz+3*cos(lambda)^3*EIz*AG-3*sin(lambda)*R^2*lambda*AE*AG-3*cos(lambda)^2*EIz*AG-3*EIz*AG*cos(lambda)+3*alpha*cos(lambda)^2*AE*EIz+5*R^2*cos(lambda)^3*AE*AG-2*R^2*cos(lambda)^4*AE*AG-3*alpha*cos(lambda)^3*AE*EIz+3*alpha*AE*EIz*cos(lambda)+3*cos(lambda)*sin(lambda)^2*EIz*AG-3*sin(lambda)*lambda*EIz*AG+5*R^2*sin(lambda)^2*cos(lambda)*AE*AG)*AE*EIz*AG/(-18*EIz^2*AG^2*lambda+36*alpha*lambda^3*AE*EIz^2*AG+36*R^2*lambda^3*AE*AG^2*EIz+36*R^2*lambda^3*AE^2*AG*alpha*EIz+18*R^4*lambda^3*AE^2*AG^2-54*R^4*lambda*AE^2*AG^2+18*alpha^2*lambda^3*AE^2*EIz^2-18*alpha^2*AE^2*EIz^2*lambda+18*lambda^3*EIz^2*AG^2-24*R^2*cos(lambda)^2*sin(lambda)*AE*AG^2*lambda^2*EIz-36*R^4*cos(lambda)^4*sin(lambda)*AE^2*AG^2-12*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2+24*cos(lambda)^5*EIz*AG^2*lambda*R^2*AE-72*cos(lambda)^2*EIz^2*AG*alpha*lambda*AE-36*cos(lambda)^4*EIz*AG^2*lambda*R^2*AE-48*cos(lambda)^3*EIz*AG^2*lambda*R^2*AE-45*cos(lambda)^2*sin(lambda)^2*EIz*AG^2*lambda*R^2*AE+36*cos(lambda)^3*sin(lambda)*EIz*AG^2*R^2*AE-108*cos(lambda)*sin(lambda)*EIz*AG^2*R^2*AE-11*R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda+36*lambda*EIz*AG^2*R^2*AE*cos(lambda)^2+96*lambda*EIz*AG^2*R^2*AE*cos(lambda)-9*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2*lambda^2-9*R^2*cos(lambda)^3*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-9*R^2*cos(lambda)^4*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz-7*R^4*cos(lambda)^4*sin(lambda)^2*AE^2*AG^2*lambda+9*R^2*cos(lambda)^4*sin(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*cos(lambda)^3*sin(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*cos(lambda)^5*sin(lambda)*AE^2*AG^2-9*R^4*sin(lambda)^3*cos(lambda)*AE^2*AG^2*lambda^2-9*R^2*sin(lambda)^3*cos(lambda)*AE^2*AG*alpha*lambda^2*EIz-24*R^2*sin(lambda)^4*AE^2*AG*lambda*alpha*cos(lambda)*EIz+24*R^4*sin(lambda)^4*AE^2*AG^2*lambda*cos(lambda)+24*R^2*sin(lambda)^4*AE*AG^2*lambda*cos(lambda)*EIz+48*R^2*cos(lambda)^3*sin(lambda)^2*AE*AG^2*lambda*EIz+45*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+36*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG-36*alpha*cos(lambda)^3*sin(lambda)*AE^2*EIz*R^2*AG+108*alpha*cos(lambda)*sin(lambda)*AE^2*EIz*R^2*AG-24*R^4*sin(lambda)^3*AE^2*AG^2*lambda^2+96*R^4*lambda*AE^2*AG^2*cos(lambda)-108*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)^2+48*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)-24*R^2*sin(lambda)^3*AE*AG^2*lambda^2*EIz-36*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)^2+6*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)+24*R^4*AE^2*AG^2*sin(lambda)*lambda^2+24*R^2*AE^2*AG*sin(lambda)*alpha*lambda^2*EIz+48*R^2*AE^2*AG*sin(lambda)^2*lambda*alpha*cos(lambda)*EIz-48*R^4*AE^2*AG^2*sin(lambda)^2*lambda*cos(lambda)-48*R^2*AE*AG^2*sin(lambda)^2*lambda*cos(lambda)*EIz+24*R^2*AE*AG^2*sin(lambda)*lambda^2*EIz+84*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)^2-78*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)-24*R^4*cos(lambda)^2*sin(lambda)*AE^2*AG^2*lambda^2-24*R^2*cos(lambda)^2*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-48*R^2*cos(lambda)^3*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz+48*R^4*cos(lambda)^3*sin(lambda)^2*AE^2*AG^2*lambda-2*R^4*AE^2*AG^2*cos(lambda)^4*lambda-48*R^4*AE^2*AG^2*cos(lambda)^3*lambda+9*R^2*lambda^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)*EIz+9*R^4*lambda^2*AE^2*AG^2*cos(lambda)*sin(lambda)+9*R^2*lambda^2*AE*AG^2*cos(lambda)*sin(lambda)*EIz+72*EIz*AG^2*R^2*AE*sin(lambda)-72*alpha*AE^2*EIz*R^2*AG*sin(lambda)-18*alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda+36*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+36*alpha^2*cos(lambda)^2*AE^2*EIz^2*lambda+48*R^2*AE^2*AG*cos(lambda)^3*lambda*alpha*EIz-8*R^4*cos(lambda)^6*AE^2*AG^2*lambda-24*R^2*cos(lambda)^5*AE^2*AG*lambda*alpha*EIz+24*R^4*cos(lambda)^5*AE^2*AG^2*lambda+R^4*sin(lambda)^4*cos(lambda)^2*AE^2*AG^2*lambda+...
   9*R^2*sin(lambda)^4*cos(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*sin(lambda)^3*cos(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*sin(lambda)^3*cos(lambda)^3*AE^2*AG^2-18*cos(lambda)^2*sin(lambda)^2*EIz^2*AG^2*lambda-8*R^4*lambda*AE^2*AG^2*cos(lambda)^2-18*alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda-24*R^2*sin(lambda)^3*AE^2*AG*alpha*lambda^2*EIz-9*R^2*sin(lambda)^4*cos(lambda)^2*AE^2*AG*lambda*alpha*EIz-36*R^2*AE*AG^2*sin(lambda)^2*lambda*EIz-36*R^4*AE^2*AG^2*sin(lambda)^2*lambda-36*R^2*AE^2*AG*sin(lambda)^2*alpha*lambda*EIz-36*R^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)^3*EIz+36*R^2*AE*AG^2*cos(lambda)*sin(lambda)^3*EIz+48*R^4*sin(lambda)^3*AE^2*AG^2+24*R^4*AE^2*AG^2*sin(lambda)-18*cos(lambda)^4*EIz^2*AG^2*lambda+36*cos(lambda)^2*EIz^2*AG^2*lambda+36*EIz^2*AG*alpha*lambda*AE-72*lambda*EIz*AG^2*R^2*AE+36*cos(lambda)^4*EIz^2*AG*lambda*alpha*AE);
k(2,1)=k(1,2);
k(2,2)=6*(6*R^2*lambda^2*AE*AG+6*alpha*lambda^2*AE*EIz-6*lambda*alpha*cos(lambda)*sin(lambda)*AE*EIz-8*lambda*R^2*sin(lambda)^3*AE*AG+8*sin(lambda)*R^2*lambda*AE*AG-8*lambda*R^2*cos(lambda)^2*sin(lambda)*AE*AG-3*lambda*R^2*cos(lambda)^3*sin(lambda)*AE*AG-3*lambda*R^2*sin(lambda)^3*cos(lambda)*AE*AG+6*lambda*cos(lambda)*sin(lambda)*EIz*AG+9*lambda*R^2*cos(lambda)*sin(lambda)*AE*AG+6*lambda^2*EIz*AG-12*R^2*AE*AG*sin(lambda)^2)*AE*EIz*AG/(R*(-18*EIz^2*AG^2*lambda+36*alpha*lambda^3*AE*EIz^2*AG+36*R^2*lambda^3*AE*AG^2*EIz+36*R^2*lambda^3*AE^2*AG*alpha*EIz+18*R^4*lambda^3*AE^2*AG^2-54*R^4*lambda*AE^2*AG^2+18*alpha^2*lambda^3*AE^2*EIz^2-18*alpha^2*AE^2*EIz^2*lambda+18*lambda^3*EIz^2*AG^2-24*R^2*cos(lambda)^2*sin(lambda)*AE*AG^2*lambda^2*EIz-36*R^4*cos(lambda)^4*sin(lambda)*AE^2*AG^2-12*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2+24*cos(lambda)^5*EIz*AG^2*lambda*R^2*AE-72*cos(lambda)^2*EIz^2*AG*alpha*lambda*AE-36*cos(lambda)^4*EIz*AG^2*lambda*R^2*AE-48*cos(lambda)^3*EIz*AG^2*lambda*R^2*AE-45*cos(lambda)^2*sin(lambda)^2*EIz*AG^2*lambda*R^2*AE+36*cos(lambda)^3*sin(lambda)*EIz*AG^2*R^2*AE-108*cos(lambda)*sin(lambda)*EIz*AG^2*R^2*AE-11*R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda+36*lambda*EIz*AG^2*R^2*AE*cos(lambda)^2+96*lambda*EIz*AG^2*R^2*AE*cos(lambda)-9*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2*lambda^2-9*R^2*cos(lambda)^3*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-9*R^2*cos(lambda)^4*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz-7*R^4*cos(lambda)^4*sin(lambda)^2*AE^2*AG^2*lambda+9*R^2*cos(lambda)^4*sin(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*cos(lambda)^3*sin(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*cos(lambda)^5*sin(lambda)*AE^2*AG^2-9*R^4*sin(lambda)^3*cos(lambda)*AE^2*AG^2*lambda^2-9*R^2*sin(lambda)^3*cos(lambda)*AE^2*AG*alpha*lambda^2*EIz-24*R^2*sin(lambda)^4*AE^2*AG*lambda*alpha*cos(lambda)*EIz+24*R^4*sin(lambda)^4*AE^2*AG^2*lambda*cos(lambda)+24*R^2*sin(lambda)^4*AE*AG^2*lambda*cos(lambda)*EIz+48*R^2*cos(lambda)^3*sin(lambda)^2*AE*AG^2*lambda*EIz+45*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+36*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG-36*alpha*cos(lambda)^3*sin(lambda)*AE^2*EIz*R^2*AG+108*alpha*cos(lambda)*sin(lambda)*AE^2*EIz*R^2*AG-24*R^4*sin(lambda)^3*AE^2*AG^2*lambda^2+96*R^4*lambda*AE^2*AG^2*cos(lambda)-108*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)^2+48*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)-24*R^2*sin(lambda)^3*AE*AG^2*lambda^2*EIz-36*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)^2+6*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)+24*R^4*AE^2*AG^2*sin(lambda)*lambda^2+24*R^2*AE^2*AG*sin(lambda)*alpha*lambda^2*EIz+48*R^2*AE^2*AG*sin(lambda)^2*lambda*alpha*cos(lambda)*EIz-48*R^4*AE^2*AG^2*sin(lambda)^2*lambda*cos(lambda)-48*R^2*AE*AG^2*sin(lambda)^2*lambda*cos(lambda)*EIz+24*R^2*AE*AG^2*sin(lambda)*lambda^2*EIz+84*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)^2-78*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)-24*R^4*cos(lambda)^2*sin(lambda)*AE^2*AG^2*lambda^2-24*R^2*cos(lambda)^2*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-48*R^2*cos(lambda)^3*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz+48*R^4*cos(lambda)^3*sin(lambda)^2*AE^2*AG^2*lambda-2*R^4*AE^2*AG^2*cos(lambda)^4*lambda-48*R^4*AE^2*AG^2*cos(lambda)^3*lambda+9*R^2*lambda^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)*EIz+9*R^4*lambda^2*AE^2*AG^2*cos(lambda)*sin(lambda)+9*R^2*lambda^2*AE*AG^2*cos(lambda)*sin(lambda)*EIz+72*EIz*AG^2*R^2*AE*sin(lambda)-72*alpha*AE^2*EIz*R^2*AG*sin(lambda)-18*alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda+36*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+36*alpha^2*cos(lambda)^2*AE^2*EIz^2*lambda+48*R^2*AE^2*AG*cos(lambda)^3*lambda*alpha*EIz-8*R^4*cos(lambda)^6*AE^2*AG^2*lambda-24*R^2*cos(lambda)^5*AE^2*AG*lambda*alpha*EIz+24*R^4*cos(lambda)^5*AE^2*AG^2*lambda+R^4*sin(lambda)^4*cos(lambda)^2*AE^2*AG^2*lambda+9*R^2*sin(lambda)^4*cos(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*sin(lambda)^3*cos(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*sin(lambda)^3*cos(lambda)^3*AE^2*AG^2-18*cos(lambda)^2*sin(lambda)^2*EIz^2*AG^2*lambda-8*R^4*lambda*AE^2*AG^2*cos(lambda)^2-18*alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda-...
   24*R^2*sin(lambda)^3*AE^2*AG*alpha*lambda^2*EIz-9*R^2*sin(lambda)^4*cos(lambda)^2*AE^2*AG*lambda*alpha*EIz-36*R^2*AE*AG^2*sin(lambda)^2*lambda*EIz-36*R^4*AE^2*AG^2*sin(lambda)^2*lambda-36*R^2*AE^2*AG*sin(lambda)^2*alpha*lambda*EIz-36*R^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)^3*EIz+36*R^2*AE*AG^2*cos(lambda)*sin(lambda)^3*EIz+48*R^4*sin(lambda)^3*AE^2*AG^2+24*R^4*AE^2*AG^2*sin(lambda)-18*cos(lambda)^4*EIz^2*AG^2*lambda+36*cos(lambda)^2*EIz^2*AG^2*lambda+36*EIz^2*AG*alpha*lambda*AE-72*lambda*EIz*AG^2*R^2*AE+36*cos(lambda)^4*EIz^2*AG*lambda*alpha*AE));
k(2,6)=6*(-6*alpha*cos(lambda)*sin(lambda)*AE*EIz+17*R^2*cos(lambda)*sin(lambda)*AE*AG+6*cos(lambda)*sin(lambda)*EIz*AG+12*R^2*lambda*AE*AG+12*lambda*EIz*AG-10*lambda*R^2*AE*AG*cos(lambda)-8*R^2*sin(lambda)^3*AE*AG+4*lambda*R^2*cos(lambda)^3*AE*AG+6*lambda*alpha*cos(lambda)^2*AE*EIz-6*lambda*R^2*AE*AG*cos(lambda)^2+4*lambda*R^2*sin(lambda)^2*cos(lambda)*AE*AG-6*lambda*cos(lambda)^2*EIz*AG-10*R^2*AE*AG*sin(lambda)-6*lambda*EIz*AG*cos(lambda)+6*sin(lambda)*alpha*AE*EIz-11*R^2*cos(lambda)^2*sin(lambda)*AE*AG+R^2*cos(lambda)^3*sin(lambda)*AE*AG+R^2*sin(lambda)^3*cos(lambda)*AE*AG-6*sin(lambda)*EIz*AG-6*lambda*alpha*AE*EIz*cos(lambda)+3*R^2*sin(lambda)^3*cos(lambda)^2*AE*AG+3*sin(lambda)*R^2*cos(lambda)^4*AE*AG)*AE*EIz*AG/(-18*EIz^2*AG^2*lambda+36*alpha*lambda^3*AE*EIz^2*AG+36*R^2*lambda^3*AE*AG^2*EIz+36*R^2*lambda^3*AE^2*AG*alpha*EIz+18*R^4*lambda^3*AE^2*AG^2-54*R^4*lambda*AE^2*AG^2+18*alpha^2*lambda^3*AE^2*EIz^2-18*alpha^2*AE^2*EIz^2*lambda+18*lambda^3*EIz^2*AG^2-24*R^2*cos(lambda)^2*sin(lambda)*AE*AG^2*lambda^2*EIz-36*R^4*cos(lambda)^4*sin(lambda)*AE^2*AG^2-12*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2+24*cos(lambda)^5*EIz*AG^2*lambda*R^2*AE-72*cos(lambda)^2*EIz^2*AG*alpha*lambda*AE-36*cos(lambda)^4*EIz*AG^2*lambda*R^2*AE-48*cos(lambda)^3*EIz*AG^2*lambda*R^2*AE-45*cos(lambda)^2*sin(lambda)^2*EIz*AG^2*lambda*R^2*AE+36*cos(lambda)^3*sin(lambda)*EIz*AG^2*R^2*AE-108*cos(lambda)*sin(lambda)*EIz*AG^2*R^2*AE-11*R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda+36*lambda*EIz*AG^2*R^2*AE*cos(lambda)^2+96*lambda*EIz*AG^2*R^2*AE*cos(lambda)-9*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2*lambda^2-9*R^2*cos(lambda)^3*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-9*R^2*cos(lambda)^4*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz-7*R^4*cos(lambda)^4*sin(lambda)^2*AE^2*AG^2*lambda+9*R^2*cos(lambda)^4*sin(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*cos(lambda)^3*sin(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*cos(lambda)^5*sin(lambda)*AE^2*AG^2-9*R^4*sin(lambda)^3*cos(lambda)*AE^2*AG^2*lambda^2-9*R^2*sin(lambda)^3*cos(lambda)*AE^2*AG*alpha*lambda^2*EIz-24*R^2*sin(lambda)^4*AE^2*AG*lambda*alpha*cos(lambda)*EIz+24*R^4*sin(lambda)^4*AE^2*AG^2*lambda*cos(lambda)+24*R^2*sin(lambda)^4*AE*AG^2*lambda*cos(lambda)*EIz+48*R^2*cos(lambda)^3*sin(lambda)^2*AE*AG^2*lambda*EIz+45*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+36*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG-36*alpha*cos(lambda)^3*sin(lambda)*AE^2*EIz*R^2*AG+108*alpha*cos(lambda)*sin(lambda)*AE^2*EIz*R^2*AG-24*R^4*sin(lambda)^3*AE^2*AG^2*lambda^2+96*R^4*lambda*AE^2*AG^2*cos(lambda)-108*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)^2+48*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)-24*R^2*sin(lambda)^3*AE*AG^2*lambda^2*EIz-36*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)^2+6*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)+24*R^4*AE^2*AG^2*sin(lambda)*lambda^2+24*R^2*AE^2*AG*sin(lambda)*alpha*lambda^2*EIz+48*R^2*AE^2*AG*sin(lambda)^2*lambda*alpha*cos(lambda)*EIz-48*R^4*AE^2*AG^2*sin(lambda)^2*lambda*cos(lambda)-48*R^2*AE*AG^2*sin(lambda)^2*lambda*cos(lambda)*EIz+24*R^2*AE*AG^2*sin(lambda)*lambda^2*EIz+84*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)^2-78*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)-24*R^4*cos(lambda)^2*sin(lambda)*AE^2*AG^2*lambda^2-24*R^2*cos(lambda)^2*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-48*R^2*cos(lambda)^3*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz+48*R^4*cos(lambda)^3*sin(lambda)^2*AE^2*AG^2*lambda-2*R^4*AE^2*AG^2*cos(lambda)^4*lambda-48*R^4*AE^2*AG^2*cos(lambda)^3*lambda+9*R^2*lambda^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)*EIz+9*R^4*lambda^2*AE^2*AG^2*cos(lambda)*sin(lambda)+9*R^2*lambda^2*AE*AG^2*cos(lambda)*sin(lambda)*EIz+72*EIz*AG^2*R^2*AE*sin(lambda)-72*alpha*AE^2*EIz*R^2*AG*sin(lambda)-18*alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda+36*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+36*alpha^2*cos(lambda)^2*AE^2*EIz^2*lambda+48*R^2*AE^2*AG*cos(lambda)^3*lambda*alpha*EIz-8*R^4*cos(lambda)^6*AE^2*AG^2*lambda-24*R^2*cos(lambda)^5*AE^2*AG*lambda*alpha*EIz+24*R^4*cos(lambda)^5*AE^2*AG^2*lambda+R^4*sin(lambda)^4*cos(lambda)^2*AE^2*AG^2*lambda+...
   9*R^2*sin(lambda)^4*cos(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*sin(lambda)^3*cos(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*sin(lambda)^3*cos(lambda)^3*AE^2*AG^2-18*cos(lambda)^2*sin(lambda)^2*EIz^2*AG^2*lambda-8*R^4*lambda*AE^2*AG^2*cos(lambda)^2-18*alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda-24*R^2*sin(lambda)^3*AE^2*AG*alpha*lambda^2*EIz-9*R^2*sin(lambda)^4*cos(lambda)^2*AE^2*AG*lambda*alpha*EIz-36*R^2*AE*AG^2*sin(lambda)^2*lambda*EIz-36*R^4*AE^2*AG^2*sin(lambda)^2*lambda-36*R^2*AE^2*AG*sin(lambda)^2*alpha*lambda*EIz-36*R^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)^3*EIz+36*R^2*AE*AG^2*cos(lambda)*sin(lambda)^3*EIz+48*R^4*sin(lambda)^3*AE^2*AG^2+24*R^4*AE^2*AG^2*sin(lambda)-18*cos(lambda)^4*EIz^2*AG^2*lambda+36*cos(lambda)^2*EIz^2*AG^2*lambda+36*EIz^2*AG*alpha*lambda*AE-72*lambda*EIz*AG^2*R^2*AE+36*cos(lambda)^4*EIz^2*AG*lambda*alpha*AE);
k(6,1)=k(1,6);
k(6,2)=k(2,6);
k(6,6)=(-9*R^2*cos(lambda)^3*sin(lambda)*AE*AG^2*lambda*EIz-9*R^2*sin(lambda)^4*cos(lambda)^2*AE^2*AG*alpha*EIz+36*cos(lambda)^2*EIz^2*AG^2-9*R^2*cos(lambda)^4*sin(lambda)^2*AE^2*AG*alpha*EIz-18*EIz^2*AG^2+36*alpha*lambda^2*AE*EIz^2*AG+72*R^2*lambda^2*AE*AG^2*EIz-36*R^2*AE^2*AG*alpha*EIz+18*alpha^2*lambda^2*AE^2*EIz^2+72*R^2*lambda^2*AE^2*AG*alpha*EIz+54*R^4*lambda^2*AE^2*AG^2+18*lambda^2*EIz^2*AG^2+36*cos(lambda)^4*EIz^2*AG*alpha*AE-18*alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2-24*R^4*sin(lambda)^3*AE^2*AG^2*lambda+24*R^4*sin(lambda)^4*AE^2*AG^2*cos(lambda)-48*R^4*AE^2*AG^2*sin(lambda)*lambda+72*R^4*AE^2*AG^2*sin(lambda)^2*cos(lambda)+48*R^4*cos(lambda)^3*sin(lambda)^2*AE^2*AG^2-7*R^4*cos(lambda)^4*sin(lambda)^2*AE^2*AG^2+R^4*sin(lambda)^4*cos(lambda)^2*AE^2*AG^2-59*R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2+24*cos(lambda)^5*EIz*AG^2*R^2*AE-72*cos(lambda)^2*EIz^2*AG*alpha*AE-48*R^2*AE*AG^2*cos(lambda)*EIz-36*cos(lambda)^4*EIz*AG^2*R^2*AE+24*cos(lambda)^3*EIz*AG^2*R^2*AE-18*cos(lambda)^4*EIz^2*AG^2+36*EIz^2*AG*alpha*AE-18*alpha^2*AE^2*EIz^2+36*R^2*AE*AG^2*EIz-18*R^4*AE^2*AG^2-68*R^4*AE^2*AG^2*cos(lambda)^2-50*R^4*AE^2*AG^2*cos(lambda)^4-27*R^4*lambda*AE^2*AG^2*cos(lambda)*sin(lambda)-24*R^2*sin(lambda)^3*AE*AG^2*lambda*EIz-18*cos(lambda)^2*sin(lambda)^2*EIz^2*AG^2+45*R^2*lambda*AE^2*AG*alpha*cos(lambda)*sin(lambda)*EIz+36*alpha^2*cos(lambda)^2*AE^2*EIz^2-8*R^4*cos(lambda)^6*AE^2*AG^2-24*R^2*cos(lambda)^2*sin(lambda)*AE^2*AG*alpha*lambda*EIz+48*R^2*cos(lambda)^3*sin(lambda)^2*AE*AG^2*EIz+36*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*AG+72*R^4*AE^2*AG^2*cos(lambda)^3+48*R^4*AE^2*AG^2*cos(lambda)+24*R^2*sin(lambda)^4*AE*AG^2*cos(lambda)*EIz-9*R^2*cos(lambda)^3*sin(lambda)*AE^2*AG*alpha*lambda*EIz-9*R^2*sin(lambda)^3*cos(lambda)*AE^2*AG*alpha*lambda*EIz-9*R^4*sin(lambda)^3*cos(lambda)*AE^2*AG^2*lambda+24*R^2*AE*AG^2*sin(lambda)^2*cos(lambda)*EIz-48*R^2*AE*AG^2*sin(lambda)*lambda*EIz-24*R^4*cos(lambda)^2*sin(lambda)*AE^2*AG^2*lambda-27*R^2*lambda*AE*AG^2*cos(lambda)*sin(lambda)*EIz+9*R^2*cos(lambda)^4*sin(lambda)^2*AE*AG^2*EIz+24*R^4*cos(lambda)^5*AE^2*AG^2-18*alpha^2*cos(lambda)^4*AE^2*EIz^2+36*alpha*cos(lambda)^4*AE^2*EIz*R^2*AG+9*R^2*sin(lambda)^4*cos(lambda)^2*AE*AG^2*EIz-45*cos(lambda)^2*sin(lambda)^2*EIz*AG^2*R^2*AE-9*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2*lambda+45*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*R^2*AG-48*R^2*AE^2*AG*sin(lambda)*alpha*lambda*EIz-24*R^2*sin(lambda)^4*AE^2*AG*alpha*cos(lambda)*EIz-24*R^2*sin(lambda)^3*AE^2*AG*alpha*lambda*EIz-24*R^2*cos(lambda)^5*AE^2*AG*alpha*EIz-9*R^2*sin(lambda)^3*cos(lambda)*AE*AG^2*lambda*EIz-48*R^2*cos(lambda)^3*sin(lambda)^2*AE^2*AG*alpha*EIz-24*R^2*cos(lambda)^2*sin(lambda)*AE*AG^2*lambda*EIz-24*R^2*AE^2*AG*sin(lambda)^2*alpha*cos(lambda)*EIz-24*R^2*AE^2*AG*cos(lambda)^3*alpha*EIz+48*R^2*AE^2*AG*cos(lambda)*alpha*EIz)*EIz/(R*(-18*EIz^2*AG^2*lambda+36*alpha*lambda^3*AE*EIz^2*AG+36*R^2*lambda^3*AE*AG^2*EIz+36*R^2*lambda^3*AE^2*AG*alpha*EIz+18*R^4*lambda^3*AE^2*AG^2-54*R^4*lambda*AE^2*AG^2+18*alpha^2*lambda^3*AE^2*EIz^2-18*alpha^2*AE^2*EIz^2*lambda+18*lambda^3*EIz^2*AG^2-24*R^2*cos(lambda)^2*sin(lambda)*AE*AG^2*lambda^2*EIz-36*R^4*cos(lambda)^4*sin(lambda)*AE^2*AG^2-12*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2+24*cos(lambda)^5*EIz*AG^2*lambda*R^2*AE-72*cos(lambda)^2*EIz^2*AG*alpha*lambda*AE-36*cos(lambda)^4*EIz*AG^2*lambda*R^2*AE-48*cos(lambda)^3*EIz*AG^2*lambda*R^2*AE-45*cos(lambda)^2*sin(lambda)^2*EIz*AG^2*lambda*R^2*AE+36*cos(lambda)^3*sin(lambda)*EIz*AG^2*R^2*AE-108*cos(lambda)*sin(lambda)*EIz*AG^2*R^2*AE-11*R^4*cos(lambda)^2*sin(lambda)^2*AE^2*AG^2*lambda+36*lambda*EIz*AG^2*R^2*AE*cos(lambda)^2+96*lambda*EIz*AG^2*R^2*AE*cos(lambda)-9*R^4*cos(lambda)^3*sin(lambda)*AE^2*AG^2*lambda^2-9*R^2*cos(lambda)^3*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-9*R^2*cos(lambda)^4*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz-7*R^4*cos(lambda)^4*sin(lambda)^2*AE^2*AG^2*lambda+9*R^2*cos(lambda)^4*sin(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*cos(lambda)^3*sin(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*cos(lambda)^5*sin(lambda)*AE^2*AG^2-...
   9*R^4*sin(lambda)^3*cos(lambda)*AE^2*AG^2*lambda^2-9*R^2*sin(lambda)^3*cos(lambda)*AE^2*AG*alpha*lambda^2*EIz-24*R^2*sin(lambda)^4*AE^2*AG*lambda*alpha*cos(lambda)*EIz+24*R^4*sin(lambda)^4*AE^2*AG^2*lambda*cos(lambda)+24*R^2*sin(lambda)^4*AE*AG^2*lambda*cos(lambda)*EIz+48*R^2*cos(lambda)^3*sin(lambda)^2*AE*AG^2*lambda*EIz+45*alpha*cos(lambda)^2*sin(lambda)^2*AE^2*EIz*lambda*R^2*AG+36*alpha*cos(lambda)^2*sin(lambda)^2*AE*EIz^2*lambda*AG-36*alpha*cos(lambda)^3*sin(lambda)*AE^2*EIz*R^2*AG+108*alpha*cos(lambda)*sin(lambda)*AE^2*EIz*R^2*AG-24*R^4*sin(lambda)^3*AE^2*AG^2*lambda^2+96*R^4*lambda*AE^2*AG^2*cos(lambda)-108*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)^2+48*alpha*lambda*AE^2*EIz*R^2*AG*cos(lambda)-24*R^2*sin(lambda)^3*AE*AG^2*lambda^2*EIz-36*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)^2+6*R^4*sin(lambda)^3*AE^2*AG^2*cos(lambda)+24*R^4*AE^2*AG^2*sin(lambda)*lambda^2+24*R^2*AE^2*AG*sin(lambda)*alpha*lambda^2*EIz+48*R^2*AE^2*AG*sin(lambda)^2*lambda*alpha*cos(lambda)*EIz-48*R^4*AE^2*AG^2*sin(lambda)^2*lambda*cos(lambda)-48*R^2*AE*AG^2*sin(lambda)^2*lambda*cos(lambda)*EIz+24*R^2*AE*AG^2*sin(lambda)*lambda^2*EIz+84*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)^2-78*R^4*AE^2*AG^2*sin(lambda)*cos(lambda)-24*R^4*cos(lambda)^2*sin(lambda)*AE^2*AG^2*lambda^2-24*R^2*cos(lambda)^2*sin(lambda)*AE^2*AG*alpha*lambda^2*EIz-48*R^2*cos(lambda)^3*sin(lambda)^2*AE^2*AG*lambda*alpha*EIz+48*R^4*cos(lambda)^3*sin(lambda)^2*AE^2*AG^2*lambda-2*R^4*AE^2*AG^2*cos(lambda)^4*lambda-48*R^4*AE^2*AG^2*cos(lambda)^3*lambda+9*R^2*lambda^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)*EIz+9*R^4*lambda^2*AE^2*AG^2*cos(lambda)*sin(lambda)+9*R^2*lambda^2*AE*AG^2*cos(lambda)*sin(lambda)*EIz+72*EIz*AG^2*R^2*AE*sin(lambda)-72*alpha*AE^2*EIz*R^2*AG*sin(lambda)-18*alpha^2*cos(lambda)^4*AE^2*EIz^2*lambda+36*alpha*cos(lambda)^4*AE^2*EIz*lambda*R^2*AG+36*alpha^2*cos(lambda)^2*AE^2*EIz^2*lambda+48*R^2*AE^2*AG*cos(lambda)^3*lambda*alpha*EIz-8*R^4*cos(lambda)^6*AE^2*AG^2*lambda-24*R^2*cos(lambda)^5*AE^2*AG*lambda*alpha*EIz+24*R^4*cos(lambda)^5*AE^2*AG^2*lambda+R^4*sin(lambda)^4*cos(lambda)^2*AE^2*AG^2*lambda+9*R^2*sin(lambda)^4*cos(lambda)^2*AE*AG^2*lambda*EIz-9*R^2*sin(lambda)^3*cos(lambda)*AE*AG^2*lambda^2*EIz+18*R^4*sin(lambda)^3*cos(lambda)^3*AE^2*AG^2-18*cos(lambda)^2*sin(lambda)^2*EIz^2*AG^2*lambda-8*R^4*lambda*AE^2*AG^2*cos(lambda)^2-18*alpha^2*cos(lambda)^2*sin(lambda)^2*AE^2*EIz^2*lambda-24*R^2*sin(lambda)^3*AE^2*AG*alpha*lambda^2*EIz-9*R^2*sin(lambda)^4*cos(lambda)^2*AE^2*AG*lambda*alpha*EIz-36*R^2*AE*AG^2*sin(lambda)^2*lambda*EIz-36*R^4*AE^2*AG^2*sin(lambda)^2*lambda-36*R^2*AE^2*AG*sin(lambda)^2*alpha*lambda*EIz-36*R^2*AE^2*AG*alpha*cos(lambda)*sin(lambda)^3*EIz+36*R^2*AE*AG^2*cos(lambda)*sin(lambda)^3*EIz+48*R^4*sin(lambda)^3*AE^2*AG^2+24*R^4*AE^2*AG^2*sin(lambda)-18*cos(lambda)^4*EIz^2*AG^2*lambda+36*cos(lambda)^2*EIz^2*AG^2*lambda+36*EIz^2*AG*alpha*lambda*AE-72*lambda*EIz*AG^2*R^2*AE+36*cos(lambda)^4*EIz^2*AG*lambda*alpha*AE));

%stiffness from the xz-plane derivation done in Maple
k(3,3)=(sin(lambda)^2*cos(lambda)^2*EIy^2-2*sin(lambda)^2*cos(lambda)^2*EIy*GJ+sin(lambda)^2*cos(lambda)^2*GJ^2-lambda^2*GJ^2-2*lambda^2*GJ*EIy-lambda^2*EIy^2+EIy^2-2*GJ*EIy-2*cos(lambda)^2*EIy^2+4*EIy*cos(lambda)^2*GJ+GJ^2-2*cos(lambda)^2*GJ^2+cos(lambda)^4*EIy^2-2*cos(lambda)^4*EIy*GJ+cos(lambda)^4*GJ^2)*AG/(R*(-2*alpha*lambda*GJ*EIy-R^2*EIy*AG+R^2*GJ*AG+R^2*lambda*GJ*AG+R^2*lambda*EIy*AG+R^2*AG*lambda^2*EIy-R^2*GJ*AG*lambda^2+R^2*lambda*AG*cos(lambda)^4*GJ+R^2*lambda*EIy*AG*cos(lambda)^4-2*R^2*lambda*AG*cos(lambda)^2*GJ-2*R^2*lambda*EIy*AG*cos(lambda)^2-R^2*GJ*AG*cos(lambda)^4*sin(lambda)^2+R^2*AG*cos(lambda)^4*EIy*sin(lambda)^2-R^2*AG*cos(lambda)^2*EIy*lambda^2+R^2*GJ*AG*cos(lambda)^2*lambda^2+R^2*GJ*AG*sin(lambda)^2*cos(lambda)^2-R^2*AG*sin(lambda)^2*cos(lambda)^2*EIy+alpha*lambda*GJ^2*sin(lambda)^2*cos(lambda)^2+alpha*lambda*EIy^2*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*cos(lambda)^4+4*alpha*lambda*GJ*EIy*cos(lambda)^2-R^2*lambda^3*AG*GJ-R^2*lambda^3*EIy*AG+alpha*lambda*GJ^2*cos(lambda)^4+alpha*lambda*EIy^2*cos(lambda)^4-2*alpha*lambda*GJ^2*cos(lambda)^2-2*alpha*lambda*EIy^2*cos(lambda)^2+3*R^2*AG*cos(lambda)^2*EIy-R^2*GJ*AG*cos(lambda)^6+alpha*lambda*GJ^2-alpha*lambda^3*GJ^2-2*alpha*lambda^3*GJ*EIy-alpha*lambda^3*EIy^2+alpha*lambda*EIy^2+R^2*lambda*AG*sin(lambda)^2*cos(lambda)^2*GJ+R^2*lambda*EIy*AG*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*sin(lambda)^2*cos(lambda)^2-3*R^2*AG*cos(lambda)^4*EIy+R^2*AG*cos(lambda)^6*EIy-3*R^2*GJ*AG*cos(lambda)^2+3*R^2*GJ*AG*cos(lambda)^4));
k(3,4)=-GJ*(cos(lambda)^4*EIy+lambda^2*EIy-sin(lambda)^2*cos(lambda)^2*GJ-2*sin(lambda)*cos(lambda)*EIy*lambda-2*cos(lambda)^2*EIy+EIy-2*lambda*EIy-cos(lambda)^4*GJ+sin(lambda)^2*cos(lambda)^2*EIy+2*cos(lambda)^2*GJ-GJ+lambda^2*GJ+2*cos(lambda)^2*EIy*lambda)*AG/(-2*alpha*lambda*GJ*EIy-R^2*EIy*AG+R^2*GJ*AG+R^2*lambda*GJ*AG+R^2*lambda*EIy*AG+R^2*AG*lambda^2*EIy-R^2*GJ*AG*lambda^2+R^2*lambda*AG*cos(lambda)^4*GJ+R^2*lambda*EIy*AG*cos(lambda)^4-2*R^2*lambda*AG*cos(lambda)^2*GJ-2*R^2*lambda*EIy*AG*cos(lambda)^2-R^2*GJ*AG*cos(lambda)^4*sin(lambda)^2+R^2*AG*cos(lambda)^4*EIy*sin(lambda)^2-R^2*AG*cos(lambda)^2*EIy*lambda^2+R^2*GJ*AG*cos(lambda)^2*lambda^2+R^2*GJ*AG*sin(lambda)^2*cos(lambda)^2-R^2*AG*sin(lambda)^2*cos(lambda)^2*EIy+alpha*lambda*GJ^2*sin(lambda)^2*cos(lambda)^2+alpha*lambda*EIy^2*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*cos(lambda)^4+4*alpha*lambda*GJ*EIy*cos(lambda)^2-R^2*lambda^3*AG*GJ-R^2*lambda^3*EIy*AG+alpha*lambda*GJ^2*cos(lambda)^4+alpha*lambda*EIy^2*cos(lambda)^4-2*alpha*lambda*GJ^2*cos(lambda)^2-2*alpha*lambda*EIy^2*cos(lambda)^2+3*R^2*AG*cos(lambda)^2*EIy-R^2*GJ*AG*cos(lambda)^6+alpha*lambda*GJ^2-alpha*lambda^3*GJ^2-2*alpha*lambda^3*GJ*EIy-alpha*lambda^3*EIy^2+alpha*lambda*EIy^2+R^2*lambda*AG*sin(lambda)^2*cos(lambda)^2*GJ+R^2*lambda*EIy*AG*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*sin(lambda)^2*cos(lambda)^2-3*R^2*AG*cos(lambda)^4*EIy+R^2*AG*cos(lambda)^6*EIy-3*R^2*GJ*AG*cos(lambda)^2+3*R^2*GJ*AG*cos(lambda)^4);
k(3,5)=-EIy*(-cos(lambda)^4*GJ-lambda^2*GJ+sin(lambda)^2*cos(lambda)^2*EIy+2*sin(lambda)*cos(lambda)*lambda*GJ+2*cos(lambda)^2*GJ-GJ-2*lambda*GJ+cos(lambda)^4*EIy-sin(lambda)^2*cos(lambda)^2*GJ-2*cos(lambda)^2*EIy+EIy-lambda^2*EIy+2*cos(lambda)^2*lambda*GJ)*AG/(-2*alpha*lambda*GJ*EIy-R^2*EIy*AG+R^2*GJ*AG+R^2*lambda*GJ*AG+R^2*lambda*EIy*AG+R^2*AG*lambda^2*EIy-R^2*GJ*AG*lambda^2+R^2*lambda*AG*cos(lambda)^4*GJ+R^2*lambda*EIy*AG*cos(lambda)^4-2*R^2*lambda*AG*cos(lambda)^2*GJ-2*R^2*lambda*EIy*AG*cos(lambda)^2-R^2*GJ*AG*cos(lambda)^4*sin(lambda)^2+R^2*AG*cos(lambda)^4*EIy*sin(lambda)^2-R^2*AG*cos(lambda)^2*EIy*lambda^2+R^2*GJ*AG*cos(lambda)^2*lambda^2+R^2*GJ*AG*sin(lambda)^2*cos(lambda)^2-R^2*AG*sin(lambda)^2*cos(lambda)^2*EIy+alpha*lambda*GJ^2*sin(lambda)^2*cos(lambda)^2+alpha*lambda*EIy^2*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*cos(lambda)^4+4*alpha*lambda*GJ*EIy*cos(lambda)^2-R^2*lambda^3*AG*GJ-R^2*lambda^3*EIy*AG+alpha*lambda*GJ^2*cos(lambda)^4+alpha*lambda*EIy^2*cos(lambda)^4-2*alpha*lambda*GJ^2*cos(lambda)^2-2*alpha*lambda*EIy^2*cos(lambda)^2+3*R^2*AG*cos(lambda)^2*EIy-R^2*GJ*AG*cos(lambda)^6+alpha*lambda*GJ^2-alpha*lambda^3*GJ^2-2*alpha*lambda^3*GJ*EIy-alpha*lambda^3*EIy^2+alpha*lambda*EIy^2+R^2*lambda*AG*sin(lambda)^2*cos(lambda)^2*GJ+R^2*lambda*EIy*AG*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*sin(lambda)^2*cos(lambda)^2-3*R^2*AG*cos(lambda)^4*EIy+R^2*AG*cos(lambda)^6*EIy-3*R^2*GJ*AG*cos(lambda)^2+3*R^2*GJ*AG*cos(lambda)^4);
k(4,3)=k(3,4);
k(4,4)=GJ*(2*R^2*lambda*EIy*AG*sin(lambda)*cos(lambda)-2*R^2*AG*lambda^2*EIy+R^2*GJ*AG*sin(lambda)^2*cos(lambda)^2+2*alpha*lambda*EIy^2*sin(lambda)*cos(lambda)-2*alpha*lambda*GJ*EIy*sin(lambda)*cos(lambda)-2*alpha*lambda^2*GJ*EIy-2*alpha*lambda^2*EIy^2-R^2*GJ*AG*lambda^2+R^2*GJ*AG-2*R^2*sin(lambda)*cos(lambda)*EIy*AG-2*R^2*GJ*AG*cos(lambda)^2+2*R^2*lambda*EIy*AG+2*R^2*AG*cos(lambda)^3*EIy*sin(lambda)+R^2*GJ*AG*cos(lambda)^4-2*R^2*lambda*EIy*AG*cos(lambda)^2)/(R*(-2*alpha*lambda*GJ*EIy-R^2*EIy*AG+R^2*GJ*AG+R^2*lambda*GJ*AG+R^2*lambda*EIy*AG+R^2*AG*lambda^2*EIy-R^2*GJ*AG*lambda^2+R^2*lambda*AG*cos(lambda)^4*GJ+R^2*lambda*EIy*AG*cos(lambda)^4-2*R^2*lambda*AG*cos(lambda)^2*GJ-2*R^2*lambda*EIy*AG*cos(lambda)^2-R^2*GJ*AG*cos(lambda)^4*sin(lambda)^2+R^2*AG*cos(lambda)^4*EIy*sin(lambda)^2-R^2*AG*cos(lambda)^2*EIy*lambda^2+R^2*GJ*AG*cos(lambda)^2*lambda^2+R^2*GJ*AG*sin(lambda)^2*cos(lambda)^2-R^2*AG*sin(lambda)^2*cos(lambda)^2*EIy+alpha*lambda*GJ^2*sin(lambda)^2*cos(lambda)^2+alpha*lambda*EIy^2*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*cos(lambda)^4+4*alpha*lambda*GJ*EIy*cos(lambda)^2-R^2*lambda^3*AG*GJ-R^2*lambda^3*EIy*AG+alpha*lambda*GJ^2*cos(lambda)^4+alpha*lambda*EIy^2*cos(lambda)^4-2*alpha*lambda*GJ^2*cos(lambda)^2-2*alpha*lambda*EIy^2*cos(lambda)^2+3*R^2*AG*cos(lambda)^2*EIy-R^2*GJ*AG*cos(lambda)^6+alpha*lambda*GJ^2-alpha*lambda^3*GJ^2-2*alpha*lambda^3*GJ*EIy-alpha*lambda^3*EIy^2+alpha*lambda*EIy^2+R^2*lambda*AG*sin(lambda)^2*cos(lambda)^2*GJ+R^2*lambda*EIy*AG*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*sin(lambda)^2*cos(lambda)^2-3*R^2*AG*cos(lambda)^4*EIy+R^2*AG*cos(lambda)^6*EIy-3*R^2*GJ*AG*cos(lambda)^2+3*R^2*GJ*AG*cos(lambda)^4));
k(4,5)=-GJ*EIy*(-R^2*sin(lambda)^2*cos(lambda)^2*AG+2*R^2*lambda*AG*sin(lambda)*cos(lambda)-R^2*lambda^2*AG+2*alpha*lambda*EIy*cos(lambda)^2-2*alpha*lambda*GJ*cos(lambda)^2+2*alpha*lambda*GJ-2*alpha*lambda*EIy+R^2*AG*cos(lambda)^4-2*R^2*AG*cos(lambda)^2+R^2*AG)/(R*(-2*alpha*lambda*GJ*EIy-R^2*EIy*AG+R^2*GJ*AG+R^2*lambda*GJ*AG+R^2*lambda*EIy*AG+R^2*AG*lambda^2*EIy-R^2*GJ*AG*lambda^2+R^2*lambda*AG*cos(lambda)^4*GJ+R^2*lambda*EIy*AG*cos(lambda)^4-2*R^2*lambda*AG*cos(lambda)^2*GJ-2*R^2*lambda*EIy*AG*cos(lambda)^2-R^2*GJ*AG*cos(lambda)^4*sin(lambda)^2+R^2*AG*cos(lambda)^4*EIy*sin(lambda)^2-R^2*AG*cos(lambda)^2*EIy*lambda^2+R^2*GJ*AG*cos(lambda)^2*lambda^2+R^2*GJ*AG*sin(lambda)^2*cos(lambda)^2-R^2*AG*sin(lambda)^2*cos(lambda)^2*EIy+alpha*lambda*GJ^2*sin(lambda)^2*cos(lambda)^2+alpha*lambda*EIy^2*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*cos(lambda)^4+4*alpha*lambda*GJ*EIy*cos(lambda)^2-R^2*lambda^3*AG*GJ-R^2*lambda^3*EIy*AG+alpha*lambda*GJ^2*cos(lambda)^4+alpha*lambda*EIy^2*cos(lambda)^4-2*alpha*lambda*GJ^2*cos(lambda)^2-2*alpha*lambda*EIy^2*cos(lambda)^2+3*R^2*AG*cos(lambda)^2*EIy-R^2*GJ*AG*cos(lambda)^6+alpha*lambda*GJ^2-alpha*lambda^3*GJ^2-2*alpha*lambda^3*GJ*EIy-alpha*lambda^3*EIy^2+alpha*lambda*EIy^2+R^2*lambda*AG*sin(lambda)^2*cos(lambda)^2*GJ+R^2*lambda*EIy*AG*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*sin(lambda)^2*cos(lambda)^2-3*R^2*AG*cos(lambda)^4*EIy+R^2*AG*cos(lambda)^6*EIy-3*R^2*GJ*AG*cos(lambda)^2+3*R^2*GJ*AG*cos(lambda)^4));
k(5,3)=k(3,5);
k(5,4)=k(4,5);
k(5,5)=-EIy*(-2*R^2*lambda*AG*sin(lambda)*cos(lambda)*GJ+2*R^2*GJ*AG*lambda^2+R^2*AG*lambda^2*EIy+2*alpha*lambda*GJ*EIy*sin(lambda)*cos(lambda)-2*alpha*lambda*GJ^2*sin(lambda)*cos(lambda)+2*alpha*lambda^2*GJ^2+2*alpha*lambda^2*GJ*EIy-R^2*AG*sin(lambda)^2*cos(lambda)^2*EIy-R^2*AG*cos(lambda)^4*EIy+2*R^2*AG*cos(lambda)^2*EIy+2*R^2*AG*cos(lambda)^3*sin(lambda)*GJ-2*R^2*lambda*AG*cos(lambda)^2*GJ-R^2*EIy*AG-2*R^2*sin(lambda)*cos(lambda)*GJ*AG+2*R^2*lambda*GJ*AG)/(R*(-2*alpha*lambda*GJ*EIy-R^2*EIy*AG+R^2*GJ*AG+R^2*lambda*GJ*AG+R^2*lambda*EIy*AG+R^2*AG*lambda^2*EIy-R^2*GJ*AG*lambda^2+R^2*lambda*AG*cos(lambda)^4*GJ+R^2*lambda*EIy*AG*cos(lambda)^4-2*R^2*lambda*AG*cos(lambda)^2*GJ-2*R^2*lambda*EIy*AG*cos(lambda)^2-R^2*GJ*AG*cos(lambda)^4*sin(lambda)^2+R^2*AG*cos(lambda)^4*EIy*sin(lambda)^2-R^2*AG*cos(lambda)^2*EIy*lambda^2+R^2*GJ*AG*cos(lambda)^2*lambda^2+R^2*GJ*AG*sin(lambda)^2*cos(lambda)^2-R^2*AG*sin(lambda)^2*cos(lambda)^2*EIy+alpha*lambda*GJ^2*sin(lambda)^2*cos(lambda)^2+alpha*lambda*EIy^2*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*cos(lambda)^4+4*alpha*lambda*GJ*EIy*cos(lambda)^2-R^2*lambda^3*AG*GJ-R^2*lambda^3*EIy*AG+alpha*lambda*GJ^2*cos(lambda)^4+alpha*lambda*EIy^2*cos(lambda)^4-2*alpha*lambda*GJ^2*cos(lambda)^2-2*alpha*lambda*EIy^2*cos(lambda)^2+3*R^2*AG*cos(lambda)^2*EIy-R^2*GJ*AG*cos(lambda)^6+alpha*lambda*GJ^2-alpha*lambda^3*GJ^2-2*alpha*lambda^3*GJ*EIy-alpha*lambda^3*EIy^2+alpha*lambda*EIy^2+R^2*lambda*AG*sin(lambda)^2*cos(lambda)^2*GJ+R^2*lambda*EIy*AG*sin(lambda)^2*cos(lambda)^2-2*alpha*lambda*GJ*EIy*sin(lambda)^2*cos(lambda)^2-3*R^2*AG*cos(lambda)^4*EIy+R^2*AG*cos(lambda)^6*EIy-3*R^2*GJ*AG*cos(lambda)^2+3*R^2*GJ*AG*cos(lambda)^4));

%Other xy-plane stiffness coefficients from equilibrium laws
%   k([1,2,6],[1,2,6])=Kxy;
   k(7,1)=-k(1,1); 
   k(8,1)=-k(2,1); %
   k(7,2)=-k(1,2); %
   k(8,2)=-k(2,2);
   k(7,6)=-k(1,6); %
   k(8,6)=-k(2,6);      
   k(7,7)=-k(7,1);    
   k(8,7)=-k(7,2); %   
   k(8,8)=-k(8,2);   
   k(8,7)=-k(7,1); %? 	
   k(12,8)=-k(12,2); 
   k(12,7)=-k(12,1); %
   
   k(12,1)=-k(1,1)*R*sin(lambda)-k(1,2)*R*(1-cos(lambda))-k(1,6);
   k(12,2)=-k(2,1)*R*sin(lambda)-k(2,2)*R*(1-cos(lambda))-k(2,6);
   k(12,6)=-k(6,1)*R*sin(lambda)-k(6,2)*R*(1-cos(lambda))-k(6,6);
   k(12,7)=-k(7,1)*R*sin(lambda)-k(7,2)*R*(1-cos(lambda))-k(7,6);
   k(12,8)=-k(8,1)*R*sin(lambda)-k(8,2)*R*(1-cos(lambda))-k(8,6);
   k(12,12)=-k(12,1)*R*sin(lambda)-k(12,2)*R*(1-cos(lambda))-k(12,6); 	
   
   k(12,1)=k(1,2)*R*sin(lambda)-k(1,1)*R*(1-cos(lambda))-k(1,6);
   k(12,2)=k(2,2)*R*sin(lambda)-k(2,1)*R*(1-cos(lambda))-k(2,6);
   k(12,6)=k(6,2)*R*sin(lambda)-k(6,1)*R*(1-cos(lambda))-k(6,6);
   k(12,7)=k(7,2)*R*sin(lambda)-k(7,1)*R*(1-cos(lambda))-k(7,6);
   k(12,8)=k(8,2)*R*sin(lambda)-k(8,1)*R*(1-cos(lambda))-k(8,6);
   k(12,12)=k(12,2)*R*sin(lambda)-k(12,1)*R*(1-cos(lambda))-k(12,6); 	
   
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
   
   k(9,3)=k(3,5)*R*sin(lambda)-k(3,4)*R*(1-cos(lambda))-k(3,3); 
   k(9,4)=k(4,5)*R*sin(lambda)-k(4,4)*R*(1-cos(lambda))-k(4,3); 
   k(9,5)=k(5,5)*R*sin(lambda)-k(5,4)*R*(1-cos(lambda))-k(5,3); 
   k(9,9)=k(9,5)*R*sin(lambda)-k(9,4)*R*(1-cos(lambda))-k(9,3);  	
   k(9,10)=k(10,5)*R*sin(lambda)-k(10,4)*R*(1-cos(lambda))-k(10,3); 
   k(9,11)=k(11,5)*R*sin(lambda)-k(11,4)*R*(1-cos(lambda))-k(11,3); 

   %fill in the symmetric part
   for i=1:12
      for j=1:12
         if k(i,j)~=0,k(j,i)=k(i,j);end
      end
   end
   
   %fill in the symmetric part
   for i=7:12
      for j=7:12
         if i==j,k(i,j)=k(i-6,j-6);
         else,   k(i,j)=-k(i-6,j-6);
         end
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

