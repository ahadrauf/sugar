% Ideal opamp
%
% Input parameters:
%   none
%
% Nodes and variables:
%   There are four nodes, all with voltage variables ('e').
%   Node 1 and 2 provide input voltages; the output voltage
%   at node 3 is A(v2-v1) with respect to the reference voltage
%   at node 4.  Since this is an ideal opamp, it has infinite
%   gain A.
function [output] = MF_opamp(flag, R, params, x, t, nodes, varargin);
switch(flag)
case 'vars'
    output.dynamic = {1,{'e'};2,{'e'};3,{'e'};4,{'e'}};
    output.branch = {'i'};
case 'K'
    output = [ 0,  0,  0,  0,  0,
               0,  0,  0,  0,  0,
               0,  0,  0,  0,  1;
               0,  0,  0,  0, -1;
               1, -1,  0,  0,  0 ];
case 'display'
    gfx_opamp2(...
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
