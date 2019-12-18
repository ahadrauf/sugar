%Nonlinear beam, stretch
function [output] = MF_beam3dnl3(flag, R, param, q, t, nodes, varargin);

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
   [F] = computeforce(flag,R,param,q);
   output = F;
   
case 'dFdx'
   [F,dFdx] = computeforce(flag,R,param,q);
   output = dFdx;
  
case 'display'

  displaybeam(q, nodes(1).pos, R, param);
  
otherwise

  output = [];

end

%============

function [F,dF] = computeforce(flag,R,param,q)

   %rotations and displacements KG=R*KL*R',FG=R*FL,qG=R*qL. 
   o = zeros(3,3); %Matrix place-holders. 
   O = zeros(6,6); %Matrix place-holders. 
   R3 = R; %Rotate initial orientation
   dR31 = rot2local(q(4),q(5),q(6))'; %Additional rotation of node1. 
   dR32 = rot2local(q(10),q(11),q(12))'; %Additional rotation of node2. 
   R = [[R3 o; o R3] O; O [R3 o; o R3]]; %12x12 rotation. 
   dR1 = [[dR31 o; o dR31] O; O [dR31 o; o dR31]]; %12x12 rotation.    
   dR2 = [[dR32 o; o dR32] O; O [dR32 o; o dR32]]; %12x12 rotation.    
   %q is size 24x1, [displacement;velocity] => [q ; qdot]
   Q = R' * q(1:12,1); %Rotated back to local.
   Q1 = R' * (dR1' * q(1:12,1)); %Rotated back to local of local due to node1.
   Q2 = R' * (dR2' * q(1:12,1)); %Rotated back to local of local due to node2.
   
   if (Q(2)-Q(8))>0
      Q=Q2;
      dR=dR2;
   else
      Q=Q1;
      dR=dR1;
   end   
   
   %geometrical and material parameters
   A = param.w * param.h; %Cross-sectional area of the beam.
   L = param.l; %Beam length, along x-axis. 
   E = param.Youngsmodulus; %Young's modulus.   
       
   %matrix elements   
   Xnl = 0.6*E*A/L^3; %stiffness in x due to y or z displacement.
   Ynl = 0.69*E*A/L^3; %stiffness in y due to y or z displacement.
   Mnl = 0.32*E*A/L^2; %moment in z or y due to y or z displacement.
   
%Xnl=0;   
%Ynl=0;
%Mnl=0;
   
   %Local forces and moments.
   F=zeros(12,1); %reset F.
   
   Y3=abs((Q(8)-Q(2))^3);
   if abs(Q(6))<abs(Q(12))      
      Mnl=Mnl;
   else
      Mnl=-Mnl;
   end
   if (abs(Q(8))>abs(Q(2)))&(Q(8)<Q(2)) %Node2 is below Node1
      F(7) = Xnl*Y3; %Fx node2
      F(8) = Ynl*Y3; %Fy node2
      F(12)= Mnl*Y3; %Mz node2
   elseif (abs(Q(8))>abs(Q(2)))&(Q(8)>Q(2)) %Node2 is above Node1
      F(7) = Xnl*Y3; %Fx node2
      F(8) =-Ynl*Y3; %Fy node2
      F(12)= Mnl*Y3; %Mz node2
   elseif (abs(Q(2))>abs(Q(8)))&(Q(2)<Q(8)) %Node1 is below Node2
      F(1) =-Xnl*Y3; %Fx node1
      F(2) = Ynl*Y3; %Fy node1
      F(6) = Mnl*Y3; %Mz node1
   elseif (abs(Q(2))>abs(Q(8)))&(Q(2)>Q(8)) %Node1 is above Node2
      F(1) =-Xnl*Y3; %Fx node1
      F(2) =-Ynl*Y3; %Fy node1
      F(6) = Mnl*Y3; %Mz node1
   end
   
%sparse(F) 

%Local forces and moments.
   dF=zeros(12,12); %reset dF.
   Y2=(Q(8)-Q(2))^2;
   if (abs(Q(8))>abs(Q(2)))&(Q(8)<Q(2)) %Node2 is below Node1
      dF(7,8)  = 3*Xnl*Y2; %Fx node2
      dF(8,8)  = 3*Ynl*Y2; %Fy node2
      dF(12,8) = 3*Mnl*Y2; %Mz node2
   elseif (abs(Q(8))>abs(Q(2)))&(Q(8)>Q(2)) %Node2 is above Node1
      dF(7,8)  = 3*Xnl*Y2; %Fx node2
      dF(8,8)  =-3*Ynl*Y2; %Fy node2
      dF(12,8) = 3*Mnl*Y2; %Mz node2
   elseif (abs(Q(2))>abs(Q(8)))&(Q(2)<Q(8)) %Node1 is below Node2
      dF(1,2)  =-3*Xnl*Y2; %Fx node1
      dF(2,2)  = 3*Ynl*Y2; %Fy node1
      dF(6,2)  = 3*Mnl*Y2; %Mz node1
   elseif (abs(Q(2))>abs(Q(8)))&(Q(2)>Q(8)) %Node1 is above Node2
      dF(1,2)  =-3*Xnl*Y2; %Fx node1
      dF(2,2)  =-3*Ynl*Y2; %Fy node1
      dF(6,2)  = 3*Mnl*Y2; %Mz node1
   end
     
   F = dR * (R * F);   
   dF = dR * (R * dF * R') * dR';   
   %a=Q(6)
   %b=Q(12)
   
   %qT = [zeros(3,1); q(4:6,1); q(7:9,1)-q(1:3,1); q(10:12,1)]; %displacement translation to origin 
%qT=q(1:12,1);   
%   qTdRR = R' * (dR' * qT); %rotation back to the local of local.

   %Nonlinear stiffness matrix, linear model   
%Xnl=0;   
%Ynl=0;
%Mnl=0;   
   Knl = [ ... 
   0     Xnl   Xnl   0     0     0,    0    -Xnl  -Xnl   0     0     0;
   0     Ynl   Ynl   0     0     0,    0    -Ynl  -Ynl   0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0    -Mnl   0     0     0,    0     0     Mnl   0     0     0;  
   0     Mnl   0     0     0     0,    0    -Mnl   0     0     0     0;   
   0    -Xnl  -Xnl   0     0     0,    0     Xnl   Xnl   0     0     0;
   0    -Ynl  -Ynl   0     0     0,    0     Ynl   Ynl   0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0    -Mnl   0     0     0,    0     0     Mnl   0     0     0;
   0     Mnl   0     0     0     0,    0    -Mnl   0     0     0     0];
   Knl = [ ... 
   0    -Xnl  -Xnl   0     0     0,    0     0     0     0     0     0;
   0     Ynl   Ynl   0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0    -Mnl   0     0     0,    0     0     0     0     0     0;  
   0     Mnl   0     0     0     0,    0     0     0     0     0     0;   
   0     0     0     0     0     0,    0     Xnl   Xnl   0     0     0;
   0     0     0     0     0     0,    0     Ynl   Ynl   0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     0     0     0     0;
   0     0     0     0     0     0,    0     0     Mnl   0     0     0;
   0     0     0     0     0     0,    0    -Mnl   0     0     0     0];
   
%F = (-1) * Knl * ([0 qTdRR(2) qTdRR(3), 0 0 0, 0 qTdRR(8) qTdRR(9), 0 0 0].^3)'; %Local of local force.
%   F = (-1) * Knl * (qTdRR.^3); %Local of local force.
%   FdRR = dR * (R * F); %Global force.   
%   output_F = FdRR;
   
%   if strcmp(flag,'F') 
%      qTdRR(7)
%   end
   
%dFdx = 3 * (-1) * Knl * diag([0 qTdRR(2) qTdRR(3), 0 0 0, 0 qTdRR(8) qTdRR(9), 0 0 0].^2); %Local of local jacobian.
%   dFdx = 3 * (-1) * Knl * diag(qTdRR.^2); %Local of local jacobian.
%   dFdxdRR = dR * (R * dFdx); %Global jacobian.   
%   output_dFdx = 0*dFdxdRR; %Global
   
 %  Xnl*3.56e-6^3
 %  Ynl*3.56e-6^3
%   Mnl*3.56e-6^3
%   
