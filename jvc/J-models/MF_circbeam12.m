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

function [output] = MF_circbeam12(flag, R, param, q, t, nodes, varargin);

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
   
   
%constants
C=cos(lambda);
S=sin(lambda);
R1C=R*(1-C);
RS=R*S;

%Make k11. 
%stiffness from the xy-plane derivation done in Maple
k(1,1)=2*AG*EIz*AE*(lambda^2*EIz*AG+R^2*lambda^2*AE*AG+alpha*lambda^2*AE*EIz+lambda*C*S*EIz*AG+lambda*R^2*C*S*AE*AG-lambda*alpha*C*S*AE*EIz-2*R^2*AE*AG+2*R^2*C^2*AE*AG)/(R*(alpha^2*lambda^3*AE^2*EIz^2+lambda^3*EIz^2*AG^2+R^4*lambda^3*AE^2*AG^2+4*R^2*C*AE*AG^2*lambda*EIz+4*R^4*C*AE^2*AG^2*lambda+4*R^2*C*AE^2*AG*alpha*lambda*EIz-2*EIz^2*AG*lambda*alpha*C^2*AE+4*EIz*AG^2*R^2*S*AE-4*EIz*AG^2*R^2*C*S*AE-2*R^2*AE^2*AG*lambda*alpha*C^2*EIz-4*alpha*AE^2*EIz*R^2*S*AG+2*EIz*AG^2*lambda*R^2*C^2*AE+EIz^2*AG^2*lambda*C^2+4*R^4*AE^2*AG^2*S-4*R^4*AE^2*AG^2*C*S+alpha^2*AE^2*EIz^2*lambda*C^2+4*alpha*AE^2*EIz*R^2*C*S*AG-6*EIz*AG^2*R^2*lambda*AE+2*lambda^3*EIz*AG^2*R^2*AE+2*lambda^3*EIz^2*AG*alpha*AE+2*R^2*lambda^3*AE^2*AG*alpha*EIz-EIz^2*AG^2*lambda+2*EIz^2*AG*alpha*lambda*AE-5*R^4*AE^2*AG^2*lambda-alpha^2*AE^2*EIz^2*lambda-2*alpha*lambda*AE^2*EIz*R^2*AG+R^4*AE^2*AG^2*lambda*C^2)); 
k(1,2)=-2*(-lambda*EIz*AG-R^2*lambda*AE*AG+alpha*lambda*AE*EIz+lambda*C^2*EIz*AG+lambda*R^2*C^2*AE*AG-lambda*alpha*C^2*AE*EIz+2*R^2*S*AE*AG-2*R^2*C*S*AE*AG)*AE*EIz*AG/(R*(alpha^2*lambda^3*AE^2*EIz^2+lambda^3*EIz^2*AG^2+R^4*lambda^3*AE^2*AG^2+4*R^2*C*AE*AG^2*lambda*EIz+4*R^4*C*AE^2*AG^2*lambda+4*R^2*C*AE^2*AG*alpha*lambda*EIz-2*EIz^2*AG*lambda*alpha*C^2*AE+4*EIz*AG^2*R^2*S*AE-4*EIz*AG^2*R^2*C*S*AE-2*R^2*AE^2*AG*lambda*alpha*C^2*EIz-4*alpha*AE^2*EIz*R^2*S*AG+2*EIz*AG^2*lambda*R^2*C^2*AE+EIz^2*AG^2*lambda*C^2+4*R^4*AE^2*AG^2*S-4*R^4*AE^2*AG^2*C*S+alpha^2*AE^2*EIz^2*lambda*C^2+4*alpha*AE^2*EIz*R^2*C*S*AG-6*EIz*AG^2*R^2*lambda*AE+2*lambda^3*EIz*AG^2*R^2*AE+2*lambda^3*EIz^2*AG*alpha*AE+2*R^2*lambda^3*AE^2*AG*alpha*EIz-EIz^2*AG^2*lambda+2*EIz^2*AG*alpha*lambda*AE-5*R^4*AE^2*AG^2*lambda-alpha^2*AE^2*EIz^2*lambda-2*alpha*lambda*AE^2*EIz*R^2*AG+R^4*AE^2*AG^2*lambda*C^2));
k(1,6)=2*(2*lambda*EIz*AG+2*R^2*lambda*AE*AG-3*R^2*S*AE*AG-S*EIz*AG-lambda*R^2*C*AE*AG-lambda*R^2*C^2*AE*AG+lambda*alpha*C^2*AE*EIz-lambda*C^2*EIz*AG+C*S*EIz*AG+3*R^2*C*S*AE*AG-alpha*C*S*AE*EIz-C*alpha*lambda*AE*EIz-C*lambda*EIz*AG+alpha*S*AE*EIz)*AE*EIz*AG/(-alpha^2*lambda^3*AE^2*EIz^2-lambda^3*EIz^2*AG^2-R^4*lambda^3*AE^2*AG^2-4*R^2*C*AE*AG^2*lambda*EIz-4*R^4*C*AE^2*AG^2*lambda-4*R^2*C*AE^2*AG*alpha*lambda*EIz+2*EIz^2*AG*lambda*alpha*C^2*AE-4*EIz*AG^2*R^2*S*AE+4*EIz*AG^2*R^2*C*S*AE+2*R^2*AE^2*AG*lambda*alpha*C^2*EIz+4*alpha*AE^2*EIz*R^2*S*AG-2*EIz*AG^2*lambda*R^2*C^2*AE-EIz^2*AG^2*lambda*C^2-4*R^4*AE^2*AG^2*S+4*R^4*AE^2*AG^2*C*S-alpha^2*AE^2*EIz^2*lambda*C^2-4*alpha*AE^2*EIz*R^2*C*S*AG+6*EIz*AG^2*R^2*lambda*AE-2*lambda^3*EIz*AG^2*R^2*AE-2*lambda^3*EIz^2*AG*alpha*AE-2*R^2*lambda^3*AE^2*AG*alpha*EIz+EIz^2*AG^2*lambda-2*EIz^2*AG*alpha*lambda*AE+5*R^4*AE^2*AG^2*lambda+alpha^2*AE^2*EIz^2*lambda+2*alpha*lambda*AE^2*EIz*R^2*AG-R^4*AE^2*AG^2*lambda*C^2);
k(2,1)=k(1,2);
k(2,2)=-2*AG*EIz*AE*(lambda*C*S*EIz*AG-lambda^2*EIz*AG-R^2*lambda^2*AE*AG-alpha*lambda^2*AE*EIz+lambda*R^2*C*S*AE*AG-lambda*alpha*C*S*AE*EIz+2*R^2*AE*AG-4*R^2*C*AE*AG+2*R^2*C^2*AE*AG)/(R*(alpha^2*lambda^3*AE^2*EIz^2+lambda^3*EIz^2*AG^2+R^4*lambda^3*AE^2*AG^2+4*R^2*C*AE*AG^2*lambda*EIz+4*R^4*C*AE^2*AG^2*lambda+4*R^2*C*AE^2*AG*alpha*lambda*EIz-2*EIz^2*AG*lambda*alpha*C^2*AE+4*EIz*AG^2*R^2*S*AE-4*EIz*AG^2*R^2*C*S*AE-2*R^2*AE^2*AG*lambda*alpha*C^2*EIz-4*alpha*AE^2*EIz*R^2*S*AG+2*EIz*AG^2*lambda*R^2*C^2*AE+EIz^2*AG^2*lambda*C^2+4*R^4*AE^2*AG^2*S-4*R^4*AE^2*AG^2*C*S+alpha^2*AE^2*EIz^2*lambda*C^2+4*alpha*AE^2*EIz*R^2*C*S*AG-6*EIz*AG^2*R^2*lambda*AE+2*lambda^3*EIz*AG^2*R^2*AE+2*lambda^3*EIz^2*AG*alpha*AE+2*R^2*lambda^3*AE^2*AG*alpha*EIz-EIz^2*AG^2*lambda+2*EIz^2*AG*alpha*lambda*AE-5*R^4*AE^2*AG^2*lambda-alpha^2*AE^2*EIz^2*lambda-2*alpha*lambda*AE^2*EIz*R^2*AG+R^4*AE^2*AG^2*lambda*C^2));
k(2,6)=-2*(-R^2*AE*AG+alpha*lambda^2*AE*EIz+R^2*lambda^2*AE*AG-alpha*AE*EIz+alpha*C^2*AE*EIz-lambda*S*EIz*AG-lambda*R^2*S*AE*AG-lambda*alpha*S*AE*EIz+EIz*AG+lambda^2*EIz*AG-3*R^2*C^2*AE*AG+4*R^2*C*AE*AG-C^2*EIz*AG+lambda*alpha*C*S*AE*EIz-lambda*C*S*EIz*AG-lambda*R^2*C*S*AE*AG)*AE*EIz*AG/(alpha^2*lambda^3*AE^2*EIz^2+lambda^3*EIz^2*AG^2+R^4*lambda^3*AE^2*AG^2+4*R^2*C*AE*AG^2*lambda*EIz+4*R^4*C*AE^2*AG^2*lambda+4*R^2*C*AE^2*AG*alpha*lambda*EIz-2*EIz^2*AG*lambda*alpha*C^2*AE+4*EIz*AG^2*R^2*S*AE-4*EIz*AG^2*R^2*C*S*AE-2*R^2*AE^2*AG*lambda*alpha*C^2*EIz-4*alpha*AE^2*EIz*R^2*S*AG+2*EIz*AG^2*lambda*R^2*C^2*AE+EIz^2*AG^2*lambda*C^2+4*R^4*AE^2*AG^2*S-4*R^4*AE^2*AG^2*C*S+alpha^2*AE^2*EIz^2*lambda*C^2+4*alpha*AE^2*EIz*R^2*C*S*AG-6*EIz*AG^2*R^2*lambda*AE+2*lambda^3*EIz*AG^2*R^2*AE+2*lambda^3*EIz^2*AG*alpha*AE+2*R^2*lambda^3*AE^2*AG*alpha*EIz-EIz^2*AG^2*lambda+2*EIz^2*AG*alpha*lambda*AE-5*R^4*AE^2*AG^2*lambda-alpha^2*AE^2*EIz^2*lambda-2*alpha*lambda*AE^2*EIz*R^2*AG+R^4*AE^2*AG^2*lambda*C^2);
k(6,1)=k(1,6);
k(6,2)=k(2,6);
k(6,6)=-(-EIz^2*AG^2*C^2-2*lambda^2*EIz^2*AG*alpha*AE-alpha^2*lambda^2*AE^2*EIz^2-2*EIz*AG^2*R^2*AE-2*EIz^2*AG*alpha*AE+EIz^2*AG^2-lambda^2*EIz^2*AG^2+R^4*AE^2*AG^2+alpha^2*AE^2*EIz^2-4*lambda^2*EIz*AG^2*R^2*AE-3*R^4*lambda^2*AE^2*AG^2+2*C*S*EIz*AG^2*R^2*lambda*AE+4*alpha*lambda*AE^2*EIz*R^2*S*AG+4*lambda*EIz*AG^2*R^2*S*AE+4*R^4*lambda*AE^2*AG^2*S+2*R^4*lambda*AE^2*AG^2*C*S-2*R^2*lambda*AE^2*AG*alpha*C*S*EIz+2*EIz^2*AG*alpha*C^2*AE+7*R^4*AE^2*AG^2*C^2-alpha^2*AE^2*EIz^2*C^2-8*R^4*AE^2*AG^2*C+2*R^2*AE^2*AG*alpha*EIz-4*R^2*lambda^2*AE^2*AG*alpha*EIz+2*C^2*EIz*AG^2*R^2*AE-2*R^2*C^2*AE^2*AG*alpha*EIz)*EIz/(R*(alpha^2*lambda^3*AE^2*EIz^2+lambda^3*EIz^2*AG^2+R^4*lambda^3*AE^2*AG^2+4*R^2*C*AE*AG^2*lambda*EIz+4*R^4*C*AE^2*AG^2*lambda+4*R^2*C*AE^2*AG*alpha*lambda*EIz-2*EIz^2*AG*lambda*alpha*C^2*AE+4*EIz*AG^2*R^2*S*AE-4*EIz*AG^2*R^2*C*S*AE-2*R^2*AE^2*AG*lambda*alpha*C^2*EIz-4*alpha*AE^2*EIz*R^2*S*AG+2*EIz*AG^2*lambda*R^2*C^2*AE+EIz^2*AG^2*lambda*C^2+4*R^4*AE^2*AG^2*S-4*R^4*AE^2*AG^2*C*S+alpha^2*AE^2*EIz^2*lambda*C^2+4*alpha*AE^2*EIz*R^2*C*S*AG-6*EIz*AG^2*R^2*lambda*AE+2*lambda^3*EIz*AG^2*R^2*AE+2*lambda^3*EIz^2*AG*alpha*AE+2*R^2*lambda^3*AE^2*AG*alpha*EIz-EIz^2*AG^2*lambda+2*EIz^2*AG*alpha*lambda*AE-5*R^4*AE^2*AG^2*lambda-alpha^2*AE^2*EIz^2*lambda-2*alpha*lambda*AE^2*EIz*R^2*AG+R^4*AE^2*AG^2*lambda*C^2));
%stiffness from the xz-plane derivation done in Maple
k(3,3)=-(-2*GJ*C^2*EIy+C^2*GJ^2+C^2*EIy^2+lambda^2*GJ^2+2*lambda^2*GJ*EIy+lambda^2*EIy^2-GJ^2+2*EIy*GJ-EIy^2)*GJ*AG/(R*(alpha*lambda*EIy^2*GJ+alpha*lambda*GJ^3-2*alpha*lambda^3*EIy*GJ^2-alpha*lambda^3*EIy^2*GJ-2*alpha*lambda*EIy*GJ^2-4*EIy*AG*R^2*C*S*GJ+2*EIy*AG*R^2*lambda*GJ*C^2+2*alpha*lambda*EIy*GJ^2*C^2-alpha*lambda*EIy^2*GJ*C^2-AG*R^2*lambda*C^2*GJ^2-EIy^2*AG*R^2*lambda*C^2-4*EIy^2*AG*R^2*S+4*EIy^2*AG*R^2*C*S-4*EIy^2*AG*R^2*C*lambda-4*EIy*AG*R^2*C*lambda*GJ+4*EIy*AG*R^2*S*GJ+2*EIy*AG*R^2*lambda*GJ+5*EIy^2*AG*R^2*lambda-EIy^2*AG*R^2*lambda^3-alpha*lambda^3*GJ^3+AG*R^2*lambda*GJ^2-AG*R^2*lambda^3*GJ^2-2*EIy*AG*R^2*lambda^3*GJ-alpha*lambda*GJ^3*C^2));
k(3,4)=2*EIy*(S*GJ+C*S*EIy+lambda*EIy-S*EIy-C*EIy*lambda+lambda*GJ-C*lambda*GJ-C*S*GJ)*GJ*AG/(-alpha*lambda*EIy^2*GJ-alpha*lambda*GJ^3+2*alpha*lambda^3*EIy*GJ^2+alpha*lambda^3*EIy^2*GJ+2*alpha*lambda*EIy*GJ^2+4*EIy*AG*R^2*C*S*GJ-2*EIy*AG*R^2*lambda*GJ*C^2-2*alpha*lambda*EIy*GJ^2*C^2+alpha*lambda*EIy^2*GJ*C^2+AG*R^2*lambda*C^2*GJ^2+EIy^2*AG*R^2*lambda*C^2+4*EIy^2*AG*R^2*S-4*EIy^2*AG*R^2*C*S+4*EIy^2*AG*R^2*C*lambda+4*EIy*AG*R^2*C*lambda*GJ-4*EIy*AG*R^2*S*GJ-2*EIy*AG*R^2*lambda*GJ-5*EIy^2*AG*R^2*lambda+EIy^2*AG*R^2*lambda^3+alpha*lambda^3*GJ^3-AG*R^2*lambda*GJ^2+AG*R^2*lambda^3*GJ^2+2*EIy*AG*R^2*lambda^3*GJ+alpha*lambda*GJ^3*C^2);
k(3,5)=-(-2*lambda^2*GJ*EIy-lambda^2*GJ^2-EIy^2-lambda^2*EIy^2+GJ^2+2*S*EIy*lambda*GJ+2*S*EIy^2*lambda-C^2*GJ^2+C^2*EIy^2)*GJ*AG/(-alpha*lambda*EIy^2*GJ-alpha*lambda*GJ^3+2*alpha*lambda^3*EIy*GJ^2+alpha*lambda^3*EIy^2*GJ+2*alpha*lambda*EIy*GJ^2+4*EIy*AG*R^2*C*S*GJ-2*EIy*AG*R^2*lambda*GJ*C^2-2*alpha*lambda*EIy*GJ^2*C^2+alpha*lambda*EIy^2*GJ*C^2+AG*R^2*lambda*C^2*GJ^2+EIy^2*AG*R^2*lambda*C^2+4*EIy^2*AG*R^2*S-4*EIy^2*AG*R^2*C*S+4*EIy^2*AG*R^2*C*lambda+4*EIy*AG*R^2*C*lambda*GJ-4*EIy*AG*R^2*S*GJ-2*EIy*AG*R^2*lambda*GJ-5*EIy^2*AG*R^2*lambda+EIy^2*AG*R^2*lambda^3+alpha*lambda^3*GJ^3-AG*R^2*lambda*GJ^2+AG*R^2*lambda^3*GJ^2+2*EIy*AG*R^2*lambda^3*GJ+alpha*lambda*GJ^3*C^2);
k(4,3)=k(3,4);
k(4,4)=2*EIy*(-AG*R^2*lambda*C*S*GJ+EIy*AG*R^2*lambda*C*S-alpha*lambda*GJ^2*C*S+alpha*lambda*EIy*GJ*C*S-2*R^2*AG*EIy+2*R^2*AG*EIy*C^2+alpha*lambda^2*GJ^2+alpha*lambda^2*EIy*GJ+AG*R^2*lambda^2*GJ+EIy*AG*R^2*lambda^2)*GJ/(R*(-alpha*lambda*EIy^2*GJ-alpha*lambda*GJ^3+2*alpha*lambda^3*EIy*GJ^2+alpha*lambda^3*EIy^2*GJ+2*alpha*lambda*EIy*GJ^2+4*EIy*AG*R^2*C*S*GJ-2*EIy*AG*R^2*lambda*GJ*C^2-2*alpha*lambda*EIy*GJ^2*C^2+alpha*lambda*EIy^2*GJ*C^2+AG*R^2*lambda*C^2*GJ^2+EIy^2*AG*R^2*lambda*C^2+4*EIy^2*AG*R^2*S-4*EIy^2*AG*R^2*C*S+4*EIy^2*AG*R^2*C*lambda+4*EIy*AG*R^2*C*lambda*GJ-4*EIy*AG*R^2*S*GJ-2*EIy*AG*R^2*lambda*GJ-5*EIy^2*AG*R^2*lambda+EIy^2*AG*R^2*lambda^3+alpha*lambda^3*GJ^3-AG*R^2*lambda*GJ^2+AG*R^2*lambda^3*GJ^2+2*EIy*AG*R^2*lambda^3*GJ+alpha*lambda*GJ^3*C^2));
k(4,5)=-2*EIy*(-alpha*lambda*EIy*GJ+alpha*lambda*GJ^2+AG*GJ*R^2*C*S-AG*R^2*lambda*GJ*C^2-alpha*lambda*GJ^2*C^2+alpha*lambda*EIy*GJ*C^2+EIy*AG*R^2*lambda*C^2+3*EIy*AG*R^2*S-3*EIy*AG*R^2*C*S+EIy*AG*R^2*C*lambda+AG*R^2*C*lambda*GJ-AG*GJ*R^2*S-2*EIy*AG*R^2*lambda)*GJ/(R*(-alpha*lambda*EIy^2*GJ-alpha*lambda*GJ^3+2*alpha*lambda^3*EIy*GJ^2+alpha*lambda^3*EIy^2*GJ+2*alpha*lambda*EIy*GJ^2+4*EIy*AG*R^2*C*S*GJ-2*EIy*AG*R^2*lambda*GJ*C^2-2*alpha*lambda*EIy*GJ^2*C^2+alpha*lambda*EIy^2*GJ*C^2+AG*R^2*lambda*C^2*GJ^2+EIy^2*AG*R^2*lambda*C^2+4*EIy^2*AG*R^2*S-4*EIy^2*AG*R^2*C*S+4*EIy^2*AG*R^2*C*lambda+4*EIy*AG*R^2*C*lambda*GJ-4*EIy*AG*R^2*S*GJ-2*EIy*AG*R^2*lambda*GJ-5*EIy^2*AG*R^2*lambda+EIy^2*AG*R^2*lambda^3+alpha*lambda^3*GJ^3-AG*R^2*lambda*GJ^2+AG*R^2*lambda^3*GJ^2+2*EIy*AG*R^2*lambda^3*GJ+alpha*lambda*GJ^3*C^2));
k(5,3)=k(3,5);
k(5,4)=k(4,5);
k(5,5)=GJ*(2*alpha*lambda*EIy^2*GJ*C*S+4*EIy*AG*R^2*S*lambda*GJ+4*EIy^2*AG*R^2*S*lambda-2*EIy*AG*R^2*lambda*C*S*GJ-2*alpha*lambda*EIy*GJ^2*C*S+2*EIy^2*AG*R^2*lambda*C*S-R^2*AG*C^2*GJ^2-AG*GJ^2*R^2*lambda^2-2*alpha*lambda^2*EIy*GJ^2-2*alpha*lambda^2*EIy^2*GJ-4*EIy*AG*R^2*lambda^2*GJ-3*EIy^2*AG*R^2*lambda^2-2*EIy*AG*R^2*GJ*C^2+2*EIy*AG*R^2*GJ+AG*R^2*GJ^2+R^2*AG*EIy^2+7*R^2*AG*C^2*EIy^2-8*R^2*AG*C*EIy^2)/(R*(alpha*lambda*EIy^2*GJ+alpha*lambda*GJ^3-2*alpha*lambda^3*EIy*GJ^2-alpha*lambda^3*EIy^2*GJ-2*alpha*lambda*EIy*GJ^2-4*EIy*AG*R^2*C*S*GJ+2*EIy*AG*R^2*lambda*GJ*C^2+2*alpha*lambda*EIy*GJ^2*C^2-alpha*lambda*EIy^2*GJ*C^2-AG*R^2*lambda*C^2*GJ^2-EIy^2*AG*R^2*lambda*C^2-4*EIy^2*AG*R^2*S+4*EIy^2*AG*R^2*C*S-4*EIy^2*AG*R^2*C*lambda-4*EIy*AG*R^2*C*lambda*GJ+4*EIy*AG*R^2*S*GJ+2*EIy*AG*R^2*lambda*GJ+5*EIy^2*AG*R^2*lambda-EIy^2*AG*R^2*lambda^3-alpha*lambda^3*GJ^3+AG*R^2*lambda*GJ^2-AG*R^2*lambda^3*GJ^2-2*EIy*AG*R^2*lambda^3*GJ-alpha*lambda*GJ^3*C^2));

%Parts of k21 and k22 that are implicitly based on k11 terms
h(4,9) =-k(3,4)-k(5,4)*R1C-k(4,4)*RS;
k(9,4) =k(4,9);
k(4,10)=-k(4,4)+k(3,4)*RS;
k(10,4) =k(4,10);
k(4,11)=-k(5,4)+k(3,4)*R1C;
k(11,4) =k(4,11);
k(5,9) =-k(3,5)-k(5,5)*R1C-k(4,5)*RS;
k(9,5) =k(5,9);
k(5,10)=-k(4,5)+k(3,5)*RS;
k(10,5) =k(5,10);
 k(3,10)=-k(4,3)+k(3,3)*RS;
 k(10,3) =k(3,10);
k(5,11)=-k(5,5)+k(3,5)*R1C;
k(11,5) =k(5,11);
 k(3,11)=-k(5,3)+k(3,3)*R1C;
 k(11,3) =k(3,11);
k(6,7) =-k(1,6)+k(6,6)*RS;
k(7,6) =k(6,7);
k(6,8) =-k(2,6)+k(6,6)*R1C;
k(8,6) =k(6,8);
k(6,12)=-k(2,6)*R1C-k(1,6)*RS-k(6,6);
k(12,6) =k(6,12);
k(10,10)=-k(4,10)+k(3,10)*RS;
k(11,10)=-k(5,10)+k(3,10)*R1C;
k(10,11)=k(11,10);
k(11,11)=-k(5,11)+k(3,11)*R1C;
k(12,1) =-k(2,1)*R1C-k(1,1)*RS-k(6,1);
k(1,12) =k(12,1);
k(12,2) =-k(2,2)*R1C-k(1,2)*RS-k(6,2);
k(2,12) =k(12,2);
h(10,9)=-h(4,9)+h(3,9)*RS;
h(9,10)=h(10,9);
h(11,9)=-h(5,9)+h(3,9)*R1C;
h(9,11)=h(11,9);

%finish k21 and k12
k(7:9,1:6)=-k(1:3,1:6);
k(1:6,7:9)=-k(7:9,1:6)';

%finish k22
k(7:9,7:9)=k(1:3,1:3);
k(7:9,10:12)=-k(1:3,10:12);
k(12,12)=k(6,6);
k(10:12,7:9)=-k(10:12,1:3);



%rotate k from vertical to horizontal
R=rot2local(0,0,-pi/2)';
o3=zeros(3);
R=[R,o3;o3,R];
o6=zeros(6);
R=[R,o6;o6,R];
k = R*k*R';

%rotate k to global coordinates
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

   displaycircbeam4(q, nodes(1).pos, R, param.w, param.h, param.radius, param.alpha);

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

