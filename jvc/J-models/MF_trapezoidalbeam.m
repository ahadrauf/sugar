% 3D small-deflection linear triangular plate model
%
% Input parameters:
%   l1, l2, angle, h, ox, oy, oz  - plate lengths, angle between ajacent sides of 1st node.
%   h - beam height (often specified in process parameters)
%   density, fluid, viscosity, Youngmodulus -
%       material parameters specified as part of the process info
%netlist syntax
%triangle p1 [a b c][l1= l2= angle= h= ox= oy= oz= ] 
%a b c    = a is 1st node, b is node of side l1, c is node of side l2.
%l1, l2   = lengths of ajacent sides to 1st node
%angle    = angle from l1 to l2
%h        = plate layer thickness
%ox oy oz = initial orientation 

function [output] = MF_trapezoidalbeam(flag, R, param, q, t, nodes, varargin);
%By: jvclark - July2001
switch(flag)
case 'vars', output.dynamic={1{'x' 'y' 'z' 'rx' 'ry' 'rz'};2{'x' 'y' 'z' 'rx' 'ry' 'rz'}};    
case 'check'
   if (~isfield(param,'density')|~isfield(param,'fluid')|~isfield(param,'viscosity')|~isfield(param,'Youngsmodulus')),output='ERR: Missing process parameters.';
   elseif (~isfield(param,'l')|~isfield(param,'wbottom')|~isfield(param,'wtop')|~isfield(param,'h')), output = 'ERR Netlist syntax: trapezoidalbeam p1 [a b c][l= wtop= wbottom= ].';
   else, output=[];
   end
case 'M', output=macromatrix('mftrapezoidalbeam2.m',param,{'A1m','A5m'},'M');
case 'D', output=macromatrix('mftrapezoidalbeam2.m',param,{'A1m','A5m'},'D');
case 'K', output=macromatrix('mftrapezoidalbeam2.m',param,{'A1m','A5m'},'K');
case 'pos', output=R*[0 param.l;0 0;0 0];
case 'F', output=[];
case 'display', displaytrapezoidalbeam(q,nodes(1).pos,param,R); 
otherwise, output=[];
end

