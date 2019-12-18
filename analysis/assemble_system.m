% Assemble mass, damping, or stiffness matrices from individual
% elements.
%
% Inputs:
%   net   - netlist structure
%   mflag - a character flagging whether the routine should
%           assemble the mass matrix ('M'), damping matrix ('D'),
%           or stiffness matrix ('K')
%   is_sp - flag whether the matrix should be generated as sparse
%           or not (default: false)
%
% Outputs:
%   M     - the assembled matrix

function [M] = assemble_system(net, mflag, is_sp);

if nargin < 3
  is_sp = 0;
end

if (is_sp)
  % The number of nonzeros is pure guesswork right now...
  M = spalloc(net.dof, net.dof, 6*net.dof); 
else
  M = zeros(net.dof);
end

for i = 1:length(net.elements)

  elt = net.elements(i);

  [Mlocal] = feval( elt.model,  mflag, elt.R, elt.parameter, [], [], [], net,i );

  if (~isempty(Mlocal))
    j = find(elt.var_ids ~= 0);
    jdx = elt.var_ids(j);
    M(jdx,jdx) = M(jdx,jdx) + Mlocal(j,j);
  end

end
