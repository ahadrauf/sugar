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

function [output] = MF_rigidlinkbeamcorner(flag, R, params, q, t, nodes, varargin);
switch(flag)
case 'vars'
  output.dynamic = {1 {'x' 'y' 'z' 'rx' 'ry' 'rz'};2 {'x' 'y' 'z' 'rx' 'ry' 'rz'}};
case 'M'   
   %geometrical params
   H = params.h;               %Beam layer thickness, along z-axis.
   W = params.w;               %Beam width, along y-axis.
   L = params.l;               %Beam length, along x-axis.
   D = params.density;         %Poly-Si density.
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
   %mass matrix 
   M11=[ ...
     a     0     0     0     0     0;
     0     b     0     0     0     n;
     0     0     c     0     m     0;
     0     0     0     d     0     0;
     0     0     m     0     e     0;     
     0     n     0     0     0     f];
   M22=[ ...
     g     0     0     0     0     0;
     0     h     0     0     0    -n;
     0     0     i     0    -m     0;
     0     0     0     j     0     0;
     0     0    -m     0     k     0;
     0    -n     0     0     0     l];
   M21=[ ...
     v     0     0     0     0     0;
     0     oo    0     0     0     p;
     0     0     q     0     r     0;
     0     0     0     s     0     0;
     0     0    -r     0     w     0;
     0     t     0     0     0     u];  
   M=D*A*L*[M11,M21';M21,M22];
   %Rotate k into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*M*R')*To'; %Return the global 12*12 stiffness matrix.
case 'D'
   %matrix params
   H = params.h;		%Beam layer thickness, along z-axis.
   W = params.w;		%Beam width, along y-axis.
   L = params.l;		%Beam length, along x-axis.
   A = H * W;		%Beam crossectional area, yz-plane.
   Mu = params.viscosity;%Viscocity of fluid environment.
   delta = params.fluid;	%Fluid layer thickness.
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
%damping matrix. This viscous model is good for planar motion.  
   D11=[ ...
     a     0     0     0     0     0;
     0     b     0     0     0     n;
     0     0     c     0     m     0;
     0     0     0     d     0     0;
     0     0     m     0     e     0;     
     0     n     0     0     0     f];
   D22=[ ...
     g     0     0     0     0     0;
     0     h     0     0     0    -n;
     0     0     i     0    -m     0;
     0     0     0     j     0     0;
     0     0    -m     0     k     0;
     0    -n     0     0     0     l];
   D21=[ ...
     v     0     0     0     0     0;
     0     oo    0     0     0     p;
     0     0     q     0     r     0;
     0     0     0     s     0     0;
     0     0    -r     0     w     0;
     0     t     0     0     0     u];  
   D=Mu*W*L/delta*[D11,D21';D21,D22];
   %Rotate k into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*D*R')*To'; %Return the global 12*12 stiffness matrix.
case 'K'
   %matrix params
   E = params.Youngsmodulus;	%Modulus of elasticity.
   H = params.h;			%Beam layer thickness, along z-axis.
   W = params.w; 		%Beam width, along y-axis.
   L = params.l;			%Beam length, along x-axis.
   A = H * W;			%Beam crossectional area, yz-plane.
   Nu = params.Poisson;		%Poisson's ratio.
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
   K11=[ ...
     a     0     0     0     0     0; 
     0     b     0     0     0     j; 
     0     0     c     0    -i     0; 
     0     0     0     d     0     0; 
     0     0    -i     0     e     0;     
     0     j     0     0     0     g];
   K22=[ ...
     a     0     0     0     0     0; 
     0     b     0     0     0    -j; 
     0     0     c     0     i     0; 
     0     0     0     d     0     0; 
     0     0     i     0     e     0; 
     0    -j     0     0     0     g];   
   K12=[ ...         
    -a     0     0     0     0     0; 
     0    -b     0     0     0     j; 
     0     0    -c     0    -i     0; 
     0     0     0    -d     0     0; 
     0     0     i     0     f     0;     
     0     -j    0     0     0     h];  
   K=[K11,K12;K12',K22];
   %Rotate k into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*K*R')*To'; %Return the global 12*12 stiffness matrix.
case 'pos'
   % Compute relative positions of beam nodes
   x=params.l; %Relative x position of node2 before rotation.
   [RL1,RL2]=rigidlinkcoordinates(params); %[x;y;z] coordinated of both rigidlinks.
   P=RL2-RL1; %difference of rigidlink coordinates.
   output=R*[0,x+P(1);0,P(2);0,P(3)]; %R*[node1postions,node2postions]
case 'display'
   [RL1,RL2,r1,r2]=rigidlinkcoordinates(params); %[x;y;z] coordinated of both rigidlinks.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations. To.
   q=To'*q; %displacement of unlinked nodes.
   Rp=(nodes(1).pos-R*RL1);     
   displaybeamcorner(q,Rp,R,params,varargin{1},varargin{2});   
case 'check'
   S=[];
   if ~isfield(params,'h'),S=strcat(S,' Semicircularbeam model needslayer thickness (h).');end
   if ~isfield(params,'w'),S=strcat(S,' Semicircularbeam model needs width (w).');end
   if ~isfield(params,'Youngsmodulus'),S=strcat(S,' Semicircularbeam model needsprocess file Youngs modulus (Youngsmodulus).');end
   if ~isfield(params,'fluid'),S=strcat(S,' Semicircularbeam model needsprocess file viscous fluid layer thickness (fluid).');end
   if ~isfield(params,'viscosity'),S=strcat(S,' Semicircularbeam model needsprocess file viscosity of fluid layer (viscosity).');end
   if ~isfield(params,'density'),S=strcat(S,' Semicircularbeam model needsprocess file material density (density).');end
   output=[S];
case 'writebin'
  fid = varargin{1};
  var_ids = varargin{2};
  if (fid < 0)
    output = 1;
  else
    writebeam(fid, nodes(1).pos, R, params.l, params.w, params.h, var_ids);
  end
otherwise
  output = [];
end

