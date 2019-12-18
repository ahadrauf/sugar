function [output] = MF_beam3d_2(flag, R, params, q, t, nodes, varargin);
switch(flag)
case 'vars'
  output.dynamic = {1 {'x' 'y' 'z' 'rx' 'ry' 'rz'}; 2 {'x' 'y' 'z' 'rx' 'ry' 'rz'}};
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
   E = params.Youngsmodulus;	%Modulus of elasticity.
   H = params.h;			%Beam layer thickness, along z-axis.
   W = params.w; 		%Beam width, along y-axis.
   L = params.l;			%Beam length, along x-axis.
   A = H * W;			%Beam crossectional area, yz-plane.

   Nu = params.Poisson;		%Poisson's ratio.
   if isfield(params,'shearmodulus')
      G=params.shearmodulus;
   else
      G=E/(2*(1+Nu));	%Shear modulus. GJ = Torsional stiffness 
   end
   
%   if (W>H)
%     J = W*H^3*(16/3-3.36*H/W*(1-(H/W)^4/12))/16;
%   else
%     J = H*W^3*(16/3-3.36*W/H*(1-(W/H)^4/12))/16;
%   end

   %Torsional stiffness (Demeter G. Fertis 1996 p78.)
   a=W/2;
   b=H/2;
   sum=0;
   for n=0:100
      kn=(2*n+1)*pi/2/a;
      sum=sum+tanh(kn*b)/(2*n+1)^5;
   end
   J=16*a^3*b*(1/3-(64*a/pi^5/b)*sum); 

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
   
   
%shear effects (J.S.Przemieniecki p79)
%   phiy=(1+12*E*Iz/G/L^2/(L*H));
%   phiz=(1+12*E*Iy/G/L^2/(L*W));  
%   b=b/phiy;
%   c=c/phiz;
%   e=e/phiz;
%   g=g/phiy;
%   j=j/phiy;
%   i=i/phiz;
%   f=f/phiz;
%   h=h/phiy;
   
   
   
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

  % Compute relative positions of beam nodes

  if (isfield(params, 'l'))
    output = R * [0 params.l;
                  0 0;
                  0 0];
  else
    output = [];
  end
  
case 'postpos'

  [params, R] = beam_postpos(params, R, [nodes(1).pos, nodes(2).pos]);
  output.params = params;
  output.R = R;

case 'F'
 
% Need dF to go with it.

   w = params.w; %width
   h = params.h; %layer thickness
%   r = rot2local(params.ox,params.oy,params.oz);
   r = R';
   o = zeros(3,3);
   R = [r o; o r]; %rotation
   A = w*h; %cross-sectional area
   sigma = params.stress; %Stress < 0 if compressive, > 0 tensile.
   gamma = params.straingradient; %Strain gradient < 0 concave down, > 0 concave up.
   alpha = params.thermalexpansion; %elemental coeffiecient of thermal expansion
   I = h^3*w/12; %Moment of inertial.
   E = params.Youngsmodulus; %Young's modulus.
   ambient = params.ambienttemperature; %Average temperature outside beam
   if isfield(params,'T')
     T = params.T; %Average beam temperature
   else
     T = ambient;
   end
   thermalstress = alpha * E * (ambient-T); %Thermalstress < 0 if compressive, > 0 tensile.
   F_axialstress = sigma * A;
   F_straingradient = E * I * gamma;
   F_thermal = thermalstress * A; 
   F_node1 = R' * [(F_axialstress+F_thermal); 0; 0; 0; F_straingradient; 0];
%   FInertial=FTranslational(m,params.accel)+FCoriolis(m,params.omega,params.vel)+FTransverse(m,params.omegadot,[nodes(1).pos;nodes(2).pos])+FCentrifugal(m,params.omega,[nodes(1).pos;nodes(2).pos]);
   
   output = [F_node1; -F_node1];
   
case 'dFdx'
   E=params.Youngsmodulus;	%Modulus of elasticity.
   H=params.h;			%Beam layer thickness, along z-axis.
   W=params.w; 		%Beam width, along y-axis.
   L=params.l;			%Beam length, along x-axis.
   A=H*W;			%Beam crossectional area, yz-plane.
   Kxx=E*A/L;
   
	k11=[0 	Kcs2z		-Kcs2y;
		  0	0			0		;
		  0	0			0		];
	k12=[0	-Kcs1z	-Kcs1y;
	     0	0			0 	  ;
	     0	0			0	  ];
	k22=[0 	Kcs1z		Kcs1y;
		  0	0			0		;
		  0	0			0		];
	J=zeros(12,12);
	J(1:3,1:3)=k11;
	J(7:9,1:3)=k12;
	J(1:3,7:9)=k12';
	J(7:9,7:9)=k22;
   output=[J];
   
   
   
   %geometry   
   L=params.l;
   w=params.w
   h=params.h;
   %transform to local coordinates
   r=[R,o3,o3,o3,o3;o3,R,o3,o3;o3,o3,R,o3;o3,o3,o3,R];
   q=r*q;
   qx1=q(1);
   qy1=q(2);
   qz1=q(3);
   qox1=q(4);
   qoy1=q(5);
   qoz1=q(6);
   qx2=q(7);
   qy2=q(8);
   qz2=q(9);
   qox2=q(10);
   qoy2=q(11);
   qoz2=q(12);   
   %Hermite polynomial coefficients
   yf0=qy1;
   yf00=qoz1;
   yfL=qy2;
   yfLL=qoz2;
   yf0L=(yfL-yf0)/L;
   yf00L=(yf0L-yf00)/L;
   yf0LL=(yfLL-yf0L)/L;
   yf00LL=(yf0LL-yf00L)/L;
   zf0=qz1;
   zf00=qoy1;
   zfL=-qz2;
   zfLL=-qoy2;
   zf0L=(zfL-zf0)/L;
   zf00L=(zf0L-zf00)/L;
   zf0LL=(zfLL-zf0L)/L;
   zf00LL=(zf0LL-zf00L)/L;
   %sample angles   
   ox=0;
   oy=0;
   oz=0;
   resolution=10;
   for k=0:resolution-1
      s1=L*(k-1)/resolution; %first sample point
      s2=L*k/resolution; %second sample point
      hx1=qx1+(1+(qx2-qx1)/L)*s1; 
      hy1=yf0+s1*(yf00+s1*(yf00L+(s1-L)*yf00LL));
      hz1=zf0+s1*(zf00+s1*(zf00L+(s1-L)*zf00LL));
      hx2=qx1+(1+(qx2-qx1)/L)*s2;
      hy2=yf0+s2*(yf00+s2*(yf00L+(s2-L)*yf00LL));
      hz2=zf0+s2*(zf00+s2*(zf00L+(s2-L)*zf00LL));
      %angle from sample point 1 to 2   
      oy=oy+atan((hz2-hz1)/(hx2-hx1));
      oz=oz+atan((hy2-hy1)/(hx2-hx1));
      ox=ox+qox1*(resolution-k)/(resolution-1)+qox2*(k-1)/(resolution-1);
   end
            
   
   
   
   
   
   
   
   
   
   
case 'display'

  displaybeam(q, nodes(1).pos, R, params.l, params.w, params.h);
  
otherwise

  output = [];

end

