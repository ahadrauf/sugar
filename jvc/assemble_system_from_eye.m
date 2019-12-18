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

function [M] = assemble_system_from_eye(net, mflag, is_sp);
if nargin<3, is_sp=0; end
if (is_sp), M = sparse(eye(net.dof)); else M = eye(net.dof); end

for i = 1:length(net.elements)
  elt = net.elements(i);
  [Mlocal] = feval( elt.model,  mflag, elt.R, elt.parameter );
  if (~isempty(Mlocal))
    j = find(elt.var_ids ~= 0);
    jdx = elt.var_ids(j);
    M(jdx,jdx) = Mlocal(j,j);
  end
end
