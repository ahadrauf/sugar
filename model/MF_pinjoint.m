% MF_pinjoint.m
% Mechanical anchor
% model function defining a pin joint, anchored only in the z-axis
function [output] = MF_pinjoint(flag, R, param, q, t, nodes, varargin);
switch(flag)
case 'vars'
    output.ground = {1 {'z'}};
case 'check'
    if (~isfield(param, 'l') | ...
        ~isfield(param, 'w') | ...
        ~isfield(param, 'h'))
        output = 'Missing length, width, or height parameters';
    else
        output = [];
    end
case 'pos'
    output = [0; 0; 0];
case 'abspos'
    if (isfield(param, 'x') & isfield(param, 'y') & isfield(param, 'z'))
        output = [param.x; param.y; param.z];
    else
        output = [];
    end
case 'display'
    q1 = zeros(12,1);
    displaybeam(q1, nodes(1).pos, R, param.l, param.w, param.h);
otherwise
    output = [];
end