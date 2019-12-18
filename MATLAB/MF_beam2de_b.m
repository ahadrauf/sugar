function [output] = MF_beam2de_b(flag, R, params, q, t, nodes, varargin)
output = [];
switch(flag)
    case 'vars'
        output.dynamic = {1 {'x' 'y' 'rz' 'e'}; 2 {'x' 'y' 'rz' 'e'}};
        output.ground = {1 {'z' 'rx' 'ry'}; 2 {'z' 'rx' 'ry'}};
    case 'check'
        if ~min(isfield(params,{'l','w','h','Rs'}))==0
            output = 'Missing parameter: l, w, h, or Rs.';
        elseif ~min(isfield(params,{'density','fluid','viscosity','Youngsmodulus','permittivity'}))==0
            output = 'Missing material properties from process file.';
        end
	case 'M'
          rho = params.density;
          l = params.l;
          A = params.w * params.h;
          % Make rotation matrix 2D
          R = R(1:2,1:2);
          R(3,3) = 1;
          M11 = R * [140,    0,     0;
                       0,  156,  22*l;
                       0, 22*l, 4*l*l] * R';
          M12 = R * [70,      0,       0;
                      0,     54,   -13*l;
                      0,   13*l,  -3*l*l] * R';
          M22 = R * [140,     0,     0;
                       0,   156, -22*l;
                       0, -22*l, 4*l*l] * R';
          a1 = rho*A*l/420;
          output = a1 * [M11  M12;
                         M12' M22];
    case 'D'
          l = params.l;
          w = params.w;
          v = params.viscosity;
          f = params.fluid; 
          a1 = v*l*w/f/420;
          % Make rotation matrix 2D
          R = R(1:2,1:2);
          R(3,3) = 1;
          D11 = R * [140,     0,     0;
                       0,   156,  22*l;
                       0,  22*l, 4*l*l] * R';
          D12 = R * [70,     0,      0;
                      0,    54,  -13*l;
                      0,  13*l, -3*l*l] * R';
          D22 = R * [140,     0,       0;
                       0,   156,   -22*l;
                       0, -22*l,   4*l*l] * R';
          output = a1 * [D11  D12;
                         D12' D22];
    case 'K'
        %mechanical beam with conductivity
            K = MF_beam2d('K', R, params);
            i = 1:3; ii = 4:6;
            K11 = K(i,i); K12 = K(i,ii); K21 = K(ii,i); K22 = K(ii,ii); Z = zeros(3,1);    
            G = params.w / (params.l * params.Rs);
            output = [ 
                K11     Z   K12     Z
                Z'      G   Z'     -G
                K21     Z   K22     Z
                Z'     -G   Z'      G];
    case 'pos' %Compute relative positions of beam nodes
          if (isfield(params, 'l'))
            output = R * [...
                0 params.l;
                0 0;
                0 0];
          else
            output = [];
          end
    case 'postpos'   
        [params, R] = beam_postpos(params, R, [nodes(1).pos, nodes(2).pos]);
        output.params = params;
        output.R = R;
    case 'display'
        %Set up width and displacement (with zero padding for z, rx, ry) and display the beam
        q1([1 2 6 7 8 12]) = [q(1:3);q(5:7)];
        %a=varargin{1},b=varargin{2}
        displaybeam(q1, nodes(1).pos, R, params.l, params.w, params.h,  varargin{1},varargin{2});
    otherwise
        output = [];
end
