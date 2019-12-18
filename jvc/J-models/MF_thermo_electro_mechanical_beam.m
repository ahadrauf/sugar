
function [output] = MF_thermo_electro_mechanical_beam(flag, R, param, q, t, nodes, element_number, varargin)
%jvclark - Sep03

%Note: Need to include conduction. For now temperature is determined by current only.

    Poisson = 0.3		%Poisson's Ratio = 0.3
    thermcond = 2.33		%Thermal conductivity Si = 2.33e-6/C
    viscosity = 1.78e-5		%Viscosity (of air) = 1,78e-5
    fluid = 2e-6		% Between the device and the substrate.
    density = 2300		%Material density = 2300 kg/m^3
    Youngsmodulus = 165e9	%Young's modulus = 1.65e11 N/m^2
    permittivity = 8.854e-12	%permittivity: C^2/(uN.um^2)=(C.s)^2/kg.um^3;
    sheetresistance = 20	%Poly-Si sheet resistance [ohm/square]
    stress = 0
    straingradient=0
    thermalexpansion=0
    ambienttemperature=0
    ox = 0
    oy = 0
    oz = 0

width = param.w;
height = param.h;
thermal_coefficient_of_resistance = param.thermal_coefficient_of_resistance;
thermal_conductivity = param.thermal_conductivity;
thermal_coefficient_of_expansion = param.thermal_coefficient_of_expansion;
sheet_resistance = param.sheet_resistance;
thermal_conductivity = param.thermal_conductivity;
electrical_conductivity = param.electrical_conductivity;
Youngs_modulus = param.Youngs_modulus;

thermal_resistance = length / thermal_conductivity / width / height;
current = V/resistance;

delta_T = current^2 * resistance(element_number) * thermal_resistance;
delta_resistance = original_resistance * thermal_coefficient_of_resistance * delta_T;
strain = thermal_coefficient_of_expansion * delta_T;
cross_sectional_area = width * height;
axial_force = cross_sectional_area * strain * Youngs_modulus;







switch(flag)
case 'vars'
    output.dynamic = {1 {'x' 'y' 'z' 'rx' 'ry' 'rz' 'e'};
                      2 {'x' 'y' 'z' 'rx' 'ry' 'rz' 'e'}};

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

case 'postpos'

  [params, R] = beam_postpos(param, R, [nodes(1).pos, nodes(2).pos]);
  output.params = params;
  output.R = R;

case 'F'
    %called by assemble_F

    element_number
    
    persistent resistance
    
    if isempty(resistance) %first time going through for this particular run
        resistance(element_number) = 4;
    else
        resistance(element_number)=resistance(element_number)+10;
    end

    resistance(element_number) 
    
% Need dF to go with it.

   w = param.w; %width
   h = param.h; %layer thickness
%   r = rot2local(param.ox,param.oy,param.oz);
   r = R';
   o = zeros(3,3);
   R = [r o; o r]; %rotation
   A = w*h; %cross-sectional area
   sigma = param.stress; %Stress < 0 if compressive, > 0 tensile.
   gamma = param.straingradient; %Strain gradient < 0 concave down, > 0 concave up.
   alpha = param.thermalexpansion; %elemental coeffiecient of thermal expansion
   I = h^3*w/12; %Moment of inertial.
   E = param.Youngsmodulus; %Young's modulus.
   ambient = param.ambienttemperature; %Average temperature outside beam
   if isfield(param,'T')
     T = param.T; %Average beam temperature
   else
     T = ambient;
   end
   thermalstress = alpha * E * (ambient-T); %Thermalstress < 0 if compressive, > 0 tensile.
   F_axialstress = sigma * A;
   F_straingradient = E * I * gamma;
   F_thermal = thermalstress * A; 
   F_node1 = R' * [(F_axialstress+F_thermal); 0; 0; 0; F_straingradient; 0];
%   FInertial=FTranslational(m,param.accel)+FCoriolis(m,param.omega,param.vel)+FTransverse(m,param.omegadot,[nodes(1).pos;nodes(2).pos])+FCentrifugal(m,param.omega,[nodes(1).pos;nodes(2).pos]);
   
   output = [F_node1; -F_node1];
   
case 'display'

  displaybeam(q, nodes(1).pos, R, param.l, param.w, param.h);
  
otherwise

  output = [];

end

