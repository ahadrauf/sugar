function [output] = MF_gap2de_c(flag, R, params, q, t, nodes, varargin)
output = [];
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
            G = params.w / (params.l * params.Rs); %Rs = resistivity
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
        closure = 0.0e-6; %Desired contact point, root axis crossing
        moc = 1e-9; %Margin of convergence due to spring-force shift
        contact_slope = -1e3 / gap; %Steep force/gap slope 
        contact_intercept = -closure * contact_slope; 
    %gap conductance
        surface_fraction = 0.01; %fractional amount of surace in contact
        surface_contact_length = 0.1e-6; %contact layer
        G = surface_fraction*params.l / (surface_contact_length * params.Rs); %Conductance
        
    switch F_dFdx
        case 'F'
            F = zeros(16,1);
            if 0
                F = node_F_1(gap,q_local,closure,moc,contact_slope,contact_intercept,electrostatic,G,F);
            end
            if 1
                F = node_F_2(gap,q_local,closure,moc,contact_slope,contact_intercept,electrostatic,G,F,params);
            end
            %Output 
                output = R * F; %Rotate from local to global  
        case 'dFdx'            
            J = zeros(16,16);                
            if 0
                J = node_dFdx_1(gap,q_local,closure,moc,contact_slope,contact_intercept,electrostatic,G,J);
            end
            if 1
                J = node_dFdx_2(gap,q_local,closure,moc,contact_slope,contact_intercept,electrostatic,G,J,params);
            end
            %Output 
                output = R * J; %Rotate from local to global                
    end        
    
%% GapDisplay
function GapDisplay(R, params, q, nodes, varargin)
    q1([1 2 6 7 8 12]) = [q(1:3);q(5:7)];
    displaybeam(q1, nodes(1).pos, R, params.l, params.w1, params.h,  0,1);
    q2([1 2 6 7 8 12]) = [q(9:11);q(13:15)];
    displaybeam(q2, nodes(3).pos,R, params.l, params.w2, params.h,  0,1);
%%    
function F = node_F_1(gap,q_local,closure,moc,contact_slope,contact_intercept,electrostatic,G,F)
    %Node A
    	GAP = (gap + q_local(2) - q_local(10));
        if GAP > closure+moc  %Gap open
        	F(2) = electrostatic*q_local(4)^2 / GAP^2; %Fy node A. q_local(4) is voltage.
            F(4) = 0; %F(4) is current
        else %Gap closed
        	F(2) = contact_slope * GAP + contact_intercept; %Fy node A
            F(4) = G * q_local(4); %F(4) is current
        end
    %Node B
    	GAP = (gap + q_local(6) - q_local(14));
        	if GAP > closure+moc  %Gap open
            	F(6) = electrostatic*q_local(8)^2 / GAP^2; %Fy node B
                F(8) = 0; %F(8) is current
            else %Gap closed
            	F(6) = contact_slope * GAP + contact_intercept; %Fy node B
                F(8) = G * q_local(8); %F(8) is current
            end
function J = node_dFdx_1(gap,q_local,closure,moc,contact_slope,contact_intercept,electrostatic,G,J)
	%Node A
    	GAP = (gap + q_local(2) - q_local(10));
        if GAP > closure+moc
        	J(2,2) = -2 * electrostatic*q_local(4)^2 / GAP^3; 
            J(2,4) = 2 * electrostatic*q_local(4) / GAP^2; 
        else
        	J(2,2) = contact_slope; %Fy on node A
            J(4,4) = -G;
        end
    %Node B
    	GAP = (gap + q_local(6) - q_local(14));
        if GAP > closure+moc
        	J(6,6) = -2 * electrostatic*q_local(8)^2 / GAP^3; 
            J(6,8) = 2 * electrostatic*q_local(8) / GAP^2; %Fy on node B
        else
        	J(6,6) = contact_slope; %Fy on node B
            J(8,8) = -G;
        end
%%
function F = node_F_2(gap,q_local,closure,moc,contact_slope,contact_intercept,electrostatic,G,F,params)
    electrostatic = -electrostatic*2;
    %F node A
        B = params.permittivity * params.l * params.h;
        gap = params.gap;
        x = q_local(2); %x displacement 
        C = 1e-12 * 1e-49;
        %F(2) = 1/2*B*q_local(4)^2/(gap+x)^2 + 1/(gap+x)^7;
        x = q_local(6); %x displacement 
        %F(6) = 1/2*B*q_local(8)^2/(gap+x)^2 + 1/(gap+x)^7;
        F(6) = -1e-6;
        GAP = gap + q_local(6);
        GAP_contact = 0.1e-6;
        if q_local(6) < -(gap-0.1e-6);
            F(6) = abs(-1e-6 + -1 + 1e-49 / (GAP)^7);
            f=F(6);
            g=GAP;
        end
            f=F(6)
            g=GAP
        %F(6) = -.001e-6 + 1/(gap+x)^7;
                
function J = node_dFdx_2(gap,q_local,closure,moc,contact_slope,contact_intercept,electrostatic,G,J,params)
        B = params.permittivity * params.l * params.h;
        gap = params.gap;
        x = q_local(2); %x displacement 
        %J(2,2) = 2 * 1/2 * B * q_local(4)^2 / (gap + x)^3 + 7/(gap + x)^8;
        x = q_local(6); %x displacement 
        %J(6,6) = 2 * 1/2 * B * q_local(8)^2 / (gap + x)^3 + 7/(gap + x)^8;
        C = 1e-12 * 1e-49;

        J(6,6) = 0;% + C*7/(gap + x)^8;
        GAP = gap + q_local(6);
        GAP_contact = 0.1e-6;
        if q_local(6) < -(gap-0.1e-6);
            F(6) = 1e-49 * 7 / (GAP)^8;
        end
        %J(6,6) = 0 + 0*7/(gap + x)^8;