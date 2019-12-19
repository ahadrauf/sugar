function [output] = MF_gap2de_c_original(flag, R, params, q, t, nodes, varargin)
output = [];
pause;
switch(flag)
    case 'vars'
        output.dynamic = {1 {'x' 'y' 'rz' 'e'} ;2 {'x' 'y' 'rz' 'e'} ;3 {'x' 'y' 'rz' 'e'} ;4 {'x' 'y' 'rz' 'e'}};
        output.ground = {1 {'z' 'rx' 'ry'}; 2 {'z' 'rx' 'ry'}; 3 {'z' 'rx' 'ry'}; 4 {'z' 'rx' 'ry'}};
    case 'check'
        if ~min(isfield(params,{'l','w1','w2','gap','h','Rs'}))==0
            output = 'Missing parameter: l, w1, w2, h, gap, or Rs.';
        elseif ~min(isfield(params,{'density','fluid','viscosity','Youngsmodulus','permittivity'}))==0
            output = 'Missing material properties from process file.';
        end
    case 'M'
        output = beam_matrices('M', R, params); 
    case 'D'
        output = beam_matrices('D', R, params);
    case 'K'
        %upper plate
            params.w = params.w1;
            K = MF_beam2d('K', R, params);
            i = 1:3; ii = 4:6; %Row and column indices
            K11 = K(i,i); K12 = K(i,ii); K21 = K(ii,i); K22 = K(ii,ii); Z = zeros(3,1);
            G = params.w / (params.l * params.Rs);
            K1 = [ 
                K11     Z   K12     Z
                Z'      G   Z'     -G
                K21     Z   K22     Z
                Z'     -G   Z'      G];
        %lower plate
            params.w = params.w2;
            K = MF_beam2d('K', R, params);
            K11 = K(i,i); K12 = K(i,ii); K21 = K(ii,i); K22 = K(ii,ii); Z = zeros(3,1);
            G = params.w / (params.l * params.Rs);
            K2 = [ 
                K11     Z   K12     Z
                Z'      G   Z'     -G
                K21     Z   K22     Z
                Z'     -G   Z'      G];
        %combined plates
            Z = zeros(8,8);
            output = [K1 Z; Z K2];    
	case 'F'
        output = compute_forces(flag, R(1:2,1:2), params, q, t);
    case 'dFdx'
        output = compute_forces(flag, R(1:2,1:2), params, q, t);
    case 'pos'
        l = params.l;
        g = params.gap + (params.w1 + params.w2)/2;
        output = R * [ ...
            0  l   0  l;
            0  0  -g -g;
            0  0   0  0];
	case 'display'
        GapDisplay(R, params, q, nodes, varargin);
end

%% beam_matrices
function [M] = beam_matrices(flag, R, params)
	if (flag == 'M') || (flag == 'D')
        display('Error: gap2de not yet defined for mass or damping matrices.');
        output = [];
    end
    
%% compute_forces
function [output] = compute_forces(F_dFdx, r, params, q_global, t)
    %Displacement
        R = eye(16); R(1:2,1:2) = r; R(5:6,5:6) = r; R(9:10,9:10) = r; R(13:14,13:14) = r; 
        q_local = R' * q_global(1:16,1); %Rotate from global to local coordinates.
    %Constants
        e0 = params.permittivity;
        gap = params.gap;
        area = params.l * params.h;         
        electrostatic = -1/2 * 1/2 * area * e0;
        closure = 0.01e-6; %Desired contact point, root axis crossing
        moc = 1e-9; %Margin of convergence due to spring-force shift
        contact_slope = -1e3 / gap; %Steep force/gap slope 
        contact_intercept = -closure * contact_slope; 
    %gap conductance
        surface_fraction = 0.01; %fractional amount of surace in contact
        surface_contact_length = 0.1e-6; %contact layer  
        G = 1/params.R_insulator; %Conductance
        F = zeros(16,1);
        J = zeros(16,16);                
        %K = 0; %0.1; %Stiffness
        gap = params.gap; %gap = 2e-6; %Initial gap between plates
        permittivity = 8.854e-12; %Permittivity of the medium
        length = 100e-6; %Overlapping length of the plates
        thickness = 20e-6; %layer thickness of the plates
        B = permittivity*length*thickness; %Constant
        V = q_local(8)^2; %V = 30; %Applied voltage difference between plates
        gap1 = (gap + q_local(2) - q_local(10) - eps);
        gap2 = (gap + q_local(6) - q_local(14) - eps);
        U_capacitor = 1/2*(-1/2*B*V^2 / closure);
        C = (U_capacitor) / closure^-6; %C = (U_spring*0 + U_capacitor) / (gap-x0).^-6;
        F(2) =     -1/2*1/2*B*V^2*(gap1)^-2   -  (6*-C*(gap1)^-7 ); %sum of forces to zero out to find x
        F(6) =     -1/2*1/2*B*V^2*(gap2)^-2   -  (6*-C*(gap2)^-7 ); %sum of forces to zero out to find x        
        J(2,2) = -1/2*2*1/2*B*V^2*(gap1)^-3 + 7*(-6*-C*(gap1)^-8 ); %dF = slope = delta_F / delta_x
        J(6,6) = -1/2*2*1/2*B*V^2*(gap2)^-3 + 7*(-6*-C*(gap2)^-8 ); %dF = slope = delta_F / delta_x        
        if gap1 <= closure*2 && gap1 >= closure*2
            G = G/10;
            F(4) = -G * q_local(4); %1e-2;%G * q_local(4); %current = G * voltage
            J(4,4) = -G;
        end
        if gap1 >= closure*2 && gap2 <= closure*2
            G = G/10;
            F(8) = -G * q_local(8); %1e-2;%G * q_local(8);
            J(8,8) = -G;
        end
        if gap1 <= closure*2 && gap1 <= closure*2
            F(4) = -G * q_local(4); %1e-2;%G * q_local(4); %current = G * voltage
            J(4,4) = -G;
            F(8) = -G * q_local(8); %1e-2;%G * q_local(8);
            J(8,8) = -G;
        end
	switch F_dFdx
        case 'F'
            output = R * F; %Rotate from local to global
%o6=output(6), gap2
        case 'dFdx'
            output = R * J; %Rotate from local to global
    end

                       
%% GapDisplay
function GapDisplay(R, params, q, nodes, varargin)
    q1([1 2 6 7 8 12]) = [q(1:3);q(5:7)];
    displaybeam(q1, nodes(1).pos, R, params.l, params.w1, params.h,  0,1);
    q2([1 2 6 7 8 12]) = [q(9:11);q(13:15)];
    displaybeam(q2, nodes(3).pos,R, params.l, params.w2, params.h,  0,1);

    
    