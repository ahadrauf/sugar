%function [output] = MF_beamNL1(flag, R, param, q, t, nodes, varargin);
% 3D stretching nonlinear beam model
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
%

function [output] = MF_beamNL1(flag, R, param, q, t, nodes, varargin);
% jvclark - April 2002 

switch(flag)
    
    case 'vars'
        output.dynamic = {1 {'x' 'y' 'z' 'rx' 'ry' 'rz'}; 2 {'x' 'y' 'z' 'rx' 'ry' 'rz'}};
      
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

       %make a local q
        R12=rotation12(R);
        q=R12*q(1:12,1);
        
        %calculate K(q) and rotate to global
        output=R12*[nonlinearK(param,q)]*R12';
      
if 0
    '////////////'
    Fnon=8*sqrt(E^2*Iz^3/(L^6*A))*sqrt(15)*lambda+68/21*sqrt(E^2*Iz^3/(L^6*A))*sqrt(15)*lambda^3+(-1/441*sqrt(E^2*Iz^3/(L^6*A))*sqrt(15))*lambda^5+(-59/339570*sqrt(E^2*Iz^3/(L^6*A))*sqrt(15))*lambda^7+127/6486480*sqrt(E^2*Iz^3/(L^6*A))*sqrt(15)*lambda^9+(-93001/66745879200*sqrt(E^2*Iz^3/(L^6*A))*sqrt(15))*lambda^11+83309/995673571200*sqrt(E^2*Iz^3/(L^6*A))*sqrt(15)*lambda^13+(-591388933/130146882108122880*sqrt(E^2*Iz^3/(L^6*A))*sqrt(15))*lambda^15
    Flin=100e-6
    Ynon=1/3*sqrt(Iz/A)*sqrt(2)*sqrt(30)*lambda+1/630*sqrt(Iz/A)*sqrt(2)*sqrt(30)*lambda^3+(-1/10584*sqrt(Iz/A)*sqrt(2)*sqrt(30))*lambda^5+131/24449040*sqrt(Iz/A)*sqrt(2)*sqrt(30)*lambda^7+(-991/3467318400*sqrt(Iz/A)*sqrt(2)*sqrt(30))*lambda^9+233/16180819200*sqrt(Iz/A)*sqrt(2)*sqrt(30)*lambda^11+(-9614219/13979256939648000*sqrt(Iz/A)*sqrt(2)*sqrt(30))*lambda^13+21893171039/702793163383863552000*sqrt(Iz/A)*sqrt(2)*sqrt(30)*lambda^15
    Ylin=Flin*(50e-6)^3/12/E/Iz
    dLnon=dL
    dLlin=sqrt(50e-6^2 + Ylin^2) - 50e-6
end

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
      % Need dF to go with it.
      w = param.w; %width
      h = param.h; %layer thickness
      %r = rot2local(param.ox,param.oy,param.oz);
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
      %FInertial=FTranslational(m,param.accel)+FCoriolis(m,param.omega,param.vel)+FTransverse(m,param.omegadot,[nodes(1).pos;nodes(2).pos])+FCentrifugal(m,param.omega,[nodes(1).pos;nodes(2).pos]);   
      output = [F_node1; -F_node1];   
      
    case 'dFdx' %This Jacobian corresponds to the nonlinear bending deflection wrt a stretched beam.               
      R12=rotation12(R);
      %stiffness matrix, linear model
        %Rotate displacement q to local coords. Label displacements.
            q=R12*q(1:12,1);
            x1=q(1); %x-displacement of node 1.
            x2=q(1+6); %x-displacement of node 2.
            y1=q(2); %y-displacement of node 1.
            y2=q(2+6); %y-displacement of node 2.
            z1=q(3); %z-displacement of node 1.
            z2=q(3+6); %z-displacement of node 2.
            rz1=q(6); %rz-displacement of node 1.
            rz2=q(6+6); %rz-displacement of node 2.
            ry1=q(5); %ry-displacement of node 1.
            ry2=q(5+6); %ry-displacement of node 2.
        %Change in beam length.
            L0=param.l; %Original beam length. It may have stretched or contracted.
            L=L0+(x2-x1); %Projected length onto the x-axis.            
            dLy=1/10*(3*rz1/(L^2)-6*y2/(L^3)+6*y1/(L^3)+3*rz2/(L^2))^2*L^5+1/4*(-6*y1/(L^2)-4*rz1/L-2*rz2/L+6*y2/(L^2))*(3*rz1/(L^2)-6*y2/(L^3)+6*y1/(L^3)+3*rz2/(L^2))*L^4+1/6*(2*rz1*(3*rz1/(L^2)-6*y2/(L^3)+6*y1/(L^3)+3*rz2/(L^2))+(-6*y1/(L^2)-4*rz1/L-2*rz2/L+6*y2/(L^2))^2)*L^3+1/2*rz1*(-6*y1/(L^2)-4*rz1/L-2*rz2/L+6*y2/(L^2))*L^2+1/2*rz1^2*L;
            dLz=1/10*(3*ry1/(L^2)-6*z2/(L^3)+6*z1/(L^3)+3*ry2/(L^2))^2*L^5+1/4*(-6*z1/(L^2)-4*ry1/L-2*ry2/L+6*z2/(L^2))*(3*ry1/(L^2)-6*z2/(L^3)+6*z1/(L^3)+3*ry2/(L^2))*L^4+1/6*(2*ry1*(3*ry1/(L^2)-6*z2/(L^3)+6*z1/(L^3)+3*ry2/(L^2))+(-6*z1/(L^2)-4*ry1/L-2*ry2/L+6*z2/(L^2))^2)*L^3+1/2*ry1*(-6*z1/(L^2)-4*ry1/L-2*ry2/L+6*z2/(L^2))*L^2+1/2*ry1^2*L;
            dL=sqrt(dLz^2+dLy^2); %Net change in length due to xy-plane and xz-plane beam deflections.
        %Axial force caused by the change in length
            E=param.Youngsmodulus; %Modulus of elastisity.
            W=param.w; %Beam width.
            H=param.h; %Layer thickness.
            A=W*H; %Cross sectional area of beam.
            Fx=E*A/L0*dL; %Axial force caused by dL.
        %Stiffness terms for the xy-plane. See derivation online.
            Iz=W^3*H/12; %Moment of inertia about z-axis.
            lambda=L0/2*sqrt(Fx/(E*Iz)); %Nonlinear coupling parameter.
            dMz_dY=0;%-(6*E*Iz/(L^2)*lambda+2/5*E*Iz/(L^2)*lambda^3+(-2/175*E*Iz/(L^2))*lambda^5+4/7875*E*Iz/(L^2)*lambda^7+(-74/3031875*E*Iz/(L^2))*lambda^9+236/197071875*E*Iz/(L^2)*lambda^11+(-11012/186232921875*E*Iz/(L^2))*lambda^13+6616/2261399765625*E*Iz/(L^2)*lambda^15+(-16772918/115794974998828125*E*Iz/(L^2))*lambda^17+112134908/15632321624841796875*E*Iz/(L^2)*lambda^19+(-58118091532/163592245803969404296875*E*Iz/(L^2))*lambda^21+1308371336/74360111729077001953125*E*Iz/(L^2)*lambda^23);
            %dFy_dY=12*E*Iz/(L^3)+24/5*E*Iz/(L^3)*lambda^2+(-4/175*E*Iz/(L^3))*lambda^4+8/7875*E*Iz/(L^3)*lambda^6+(-148/3031875*E*Iz/(L^3))*lambda^8+472/197071875*E*Iz/(L^3)*lambda^10+(-22024/186232921875*E*Iz/(L^3))*lambda^12+13232/2261399765625*E*Iz/(L^3)*lambda^14+(-33545836/115794974998828125*E*Iz/(L^3))*lambda^16+224269816/15632321624841796875*E*Iz/(L^3)*lambda^18+(-116236183064/163592245803969404296875*E*Iz/(L^3))*lambda^20+2616742672/74360111729077001953125*E*Iz/(L^3)*lambda^22;
            dFy_dY=-(24/5*E*Iz/(L^3)*lambda^2+(-4/175*E*Iz/(L^3))*lambda^4+8/7875*E*Iz/(L^3)*lambda^6+(-148/3031875*E*Iz/(L^3))*lambda^8+472/197071875*E*Iz/(L^3)*lambda^10+(-22024/186232921875*E*Iz/(L^3))*lambda^12+13232/2261399765625*E*Iz/(L^3)*lambda^14+(-33545836/115794974998828125*E*Iz/(L^3))*lambda^16+224269816/15632321624841796875*E*Iz/(L^3)*lambda^18+(-116236183064/163592245803969404296875*E*Iz/(L^3))*lambda^20+2616742672/74360111729077001953125*E*Iz/(L^3)*lambda^22);
        %Stiffness terms for the xz-plane. See derivation online.
            Iy=H^3*W/12; %Moment of inertia about z-axis.
            lambda=L0/2*sqrt(Fx/(E*Iy)); %Nonlinear coupling parameter.
            %dMy_dZ=2*E*Iy*lambda^3*tanh(lambda)/(L^2*(lambda-tanh(lambda)));
            dFz_dZ=-(24/5*E*Iy/(L^3)*lambda^2+(-4/175*E*Iy/(L^3))*lambda^4+8/7875*E*Iy/(L^3)*lambda^6+(-148/3031875*E*Iy/(L^3))*lambda^8+472/197071875*E*Iy/(L^3)*lambda^10+(-22024/186232921875*E*Iy/(L^3))*lambda^12+13232/2261399765625*E*Iy/(L^3)*lambda^14+(-33545836/115794974998828125*E*Iy/(L^3))*lambda^16+224269816/15632321624841796875*E*Iy/(L^3)*lambda^18+(-116236183064/163592245803969404296875*E*Iy/(L^3))*lambda^20+2616742672/74360111729077001953125*E*Iy/(L^3)*lambda^22);
        %Nonlinear corrections
            J=zeros(12);
            J(2,2)    =  dFy_dY;
            J(2+6,2+6)=  dFy_dY;
            J(2,2+6)  = -dFy_dY;
            J(2+6,2)  = -dFy_dY;

            J(3,3)    =  dFz_dZ;
            J(3+6,3+6)=  dFz_dZ;            
            J(3,3+6)  = -dFz_dZ;
            J(3+6,3)  = -dFz_dZ;

        output = R12*J*R12';        
   
    case 'dKqdx'

        %make a local q
        R12=rotation12(R);
        q0=R12*q(1:12,1);
        
        %Approximate dKqdx
        d=7; %half of the number of significant floating point digits
        e=10^(-d); %small change parameter
        dKqdx=zeros(12,12); %init matrix
        for i=1:12 %for each column
            q1=q0;q2=q0; %re-init
            q1(i)=q0(i)-e;
            q2(i)=q0(i)+e;
            Kq1 = nonlinearK(param,q1) * q1(1:12,1);
            Kq2 = nonlinearK(param,q2) * q2(1:12,1);
            dKqdx(1:12,i)=(Kq2-Kq1)/(2*e);
        end      
        output = R12*[dKqdx]*R12';

    case 'display'
        displaybeam(q, nodes(1).pos, R, param.l, param.w, param.h);  
     
    otherwise
        output = [];
     
end %switch
  
%=================================================================================

function K=nonlinearK(param,q)
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
        %displacements
            x1=q(1); %x-displacement of node 1.
            x2=q(1+6); %x-displacement of node 2.
            y1=q(2); %y-displacement of node 1.
            y2=q(2+6); %y-displacement of node 2.
            z1=q(3); %z-displacement of node 1.
            z2=q(3+6); %z-displacement of node 2.
            rz1=q(6); %rz-displacement of node 1.
            rz2=q(6+6); %rz-displacement of node 2.
            ry1=q(5); %ry-displacement of node 1.
            ry2=q(5+6); %ry-displacement of node 2.
        %Change in beam length.
            L0=param.l; %Original beam length. It may have stretched or contracted.
            L=L0+(x2-x1); %Projected length onto the x-axis.            
            dLy=1/10*(3*rz1/(L^2)-6*y2/(L^3)+6*y1/(L^3)+3*rz2/(L^2))^2*L^5+1/4*(-6*y1/(L^2)-4*rz1/L-2*rz2/L+6*y2/(L^2))*(3*rz1/(L^2)-6*y2/(L^3)+6*y1/(L^3)+3*rz2/(L^2))*L^4+1/6*(2*rz1*(3*rz1/(L^2)-6*y2/(L^3)+6*y1/(L^3)+3*rz2/(L^2))+(-6*y1/(L^2)-4*rz1/L-2*rz2/L+6*y2/(L^2))^2)*L^3+1/2*rz1*(-6*y1/(L^2)-4*rz1/L-2*rz2/L+6*y2/(L^2))*L^2+1/2*rz1^2*L;
            dLz=1/10*(3*ry1/(L^2)-6*z2/(L^3)+6*z1/(L^3)+3*ry2/(L^2))^2*L^5+1/4*(-6*z1/(L^2)-4*ry1/L-2*ry2/L+6*z2/(L^2))*(3*ry1/(L^2)-6*z2/(L^3)+6*z1/(L^3)+3*ry2/(L^2))*L^4+1/6*(2*ry1*(3*ry1/(L^2)-6*z2/(L^3)+6*z1/(L^3)+3*ry2/(L^2))+(-6*z1/(L^2)-4*ry1/L-2*ry2/L+6*z2/(L^2))^2)*L^3+1/2*ry1*(-6*z1/(L^2)-4*ry1/L-2*ry2/L+6*z2/(L^2))*L^2+1/2*ry1^2*L;
            dL=sqrt(dLz^2+dLy^2); %Net change in length due to xy-plane and xz-plane beam deflections.
        %Axial force caused by the change in length
            A=W*H; %Cross sectional area of beam.
            Fx=E*A/L0*dL; %Axial force caused by dL.
        %Stiffness terms for the xy-plane. See derivation online.
            Iz=W^3*H/12; %Moment of inertia about z-axis.
            lambda=L0/2*sqrt(Fx/(E*Iz)); %Nonlinear coupling parameter.
            %dMz_dY=6*E*Iz/(L^2)*lambda+2/5*E*Iz/(L^2)*lambda^3+(-2/175*E*Iz/(L^2))*lambda^5+4/7875*E*Iz/(L^2)*lambda^7+(-74/3031875*E*Iz/(L^2))*lambda^9+236/197071875*E*Iz/(L^2)*lambda^11+(-11012/186232921875*E*Iz/(L^2))*lambda^13+6616/2261399765625*E*Iz/(L^2)*lambda^15+(-16772918/115794974998828125*E*Iz/(L^2))*lambda^17+112134908/15632321624841796875*E*Iz/(L^2)*lambda^19+(-58118091532/163592245803969404296875*E*Iz/(L^2))*lambda^21+1308371336/74360111729077001953125*E*Iz/(L^2)*lambda^23;
            dFy_dY=12*E*Iz/(L^3)+24/5*E*Iz/(L^3)*lambda^2+(-4/175*E*Iz/(L^3))*lambda^4+8/7875*E*Iz/(L^3)*lambda^6+(-148/3031875*E*Iz/(L^3))*lambda^8+472/197071875*E*Iz/(L^3)*lambda^10+(-22024/186232921875*E*Iz/(L^3))*lambda^12+13232/2261399765625*E*Iz/(L^3)*lambda^14+(-33545836/115794974998828125*E*Iz/(L^3))*lambda^16+224269816/15632321624841796875*E*Iz/(L^3)*lambda^18+(-116236183064/163592245803969404296875*E*Iz/(L^3))*lambda^20+2616742672/74360111729077001953125*E*Iz/(L^3)*lambda^22;
        %Stiffness terms for the xz-plane. See derivation online.
            Iy=H^3*W/12; %Moment of inertia about z-axis.
            lambda=L0/2*sqrt(Fx/(E*Iy)); %Nonlinear coupling parameter.
            %dMy_dZ=2*E*Iy*lambda^3*tanh(lambda)/(L^2*(lambda-tanh(lambda)));
            dFz_dZ=12*E*Iy/(L^3)+24/5*E*Iy/(L^3)*lambda^2+(-4/175*E*Iy/(L^3))*lambda^4+8/7875*E*Iy/(L^3)*lambda^6+(-148/3031875*E*Iy/(L^3))*lambda^8+472/197071875*E*Iy/(L^3)*lambda^10+(-22024/186232921875*E*Iy/(L^3))*lambda^12+13232/2261399765625*E*Iy/(L^3)*lambda^14+(-33545836/115794974998828125*E*Iy/(L^3))*lambda^16+224269816/15632321624841796875*E*Iy/(L^3)*lambda^18+(-116236183064/163592245803969404296875*E*Iy/(L^3))*lambda^20+2616742672/74360111729077001953125*E*Iy/(L^3)*lambda^22;

        %Nonlinear corrections
            b=dFy_dY;
            c=dFz_dZ;
            K11 = [ ...
                a     0     0     0     0     0; 
                0     b     0     0     0     j; 
                0     0     c     0    -i     0; 
                0     0     0     d     0     0; 
                0     0    -i     0     e     0;     
                0     j     0     0     0     g];
            K22 = [ ...
                a     0     0     0     0     0; 
                0     b     0     0     0    -j; 
                0     0     c     0     i     0; 
                0     0     0     d     0     0; 
                0     0     i     0     e     0; 
                0    -j     0     0     0     g];   
            K12 = [ ...         
               -a     0     0     0     0     0; 
                0    -b     0     0     0     j; 
                0     0    -c     0    -i     0; 
                0     0     0    -d     0     0; 
                0     0     i     0     f     0;     
                0     -j    0     0     0     h];          
        K = [K11  K12; K12' K22];
%=================================================================================

