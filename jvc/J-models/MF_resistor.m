function [output]=MF_resistor(flag, R, params, q, t, nodes, varargin);
switch(flag)
case 'vars'
  output.dynamic = {1 {'e'}; 2 {'e'}};
case 'check'
    if (~isfield(params, 'R') & ~isfield(params, 'G'))
        output = 'Unspecified resistance/conductance';
    else
        output = [];
    end
case 'K'
    if (isfield(params, 'R'))
        G = 1/params.R;
    else
        G = params.G;
    end
    output = [G,-G;-G,G];    
case 'pos' %Compute relative positions of nodes 
    if (isfield(params, 'len'))
        v = (rot2global(params.ox1,params.oy1,params.oz1)*[params.L1;0;0] + ...
            [params.len;0;0] + ...
            rot2global(params.ox2,params.oy2,params.oz2)*[params.L2;0;0]);
        output = R*[0,v(1);0,v(2);0,v(3)];
    else
        output = [];
    end
case 'postpos'
    [params, R] = postpos2(params, R, [nodes(1).pos, nodes(2).pos]);
    output.params = params;
    output.R = R;
case 'display'
    gfx_resistor(...
        nodes(1).pos(1), ...
        nodes(1).pos(2), ...
        nodes(1).pos(3), ...
        R, ...
        rot2global(params.ox1,params.oy1,params.oz1), ...
        rot2global(params.ox2,params.oy2,params.oz2), ...
        params.len, ...
        params.L1, ...
        params.L2, ...
        params.linewidth, ...
        params.nodewidth);
otherwise
    output = [];
end

