% Assemble the Nonlinear K, a variable stiffness matrix
%
% Inputs:
%   net   - netlist structure
%   q     - state vector of [x;xdot]
%           (note: xdot can be ommitted, in which case it is presumed 0)
%   t     - time
%   is_sp - (optional) flag whether to return Jacobian in sparse format.
%           default is false.
%
% Output:
%   NonlinearK - a variable stiffness matrix.

function [NonlinearK] = assemble_NonlinearK(net, q, t, is_sp);
if (nargin < 4)
  is_sp = 0;
end

if (is_sp)
  % Go to sparse form... pure guess work on nnz
  NonlinearK = spalloc(net.dof, net.dof, 6*net.dof);
else
  NonlinearK = zeros(net.dof);
end

has_xdot = (length(q) == 2*net.dof);

for i = 1:length(net.elements)
  elt = net.elements(i);
  if ~isempty(elt.var_ids)

    % Find the local indices of the ungrounded variables (j), 
    % and the corresponding global indices (jdx)

    j = find(elt.var_ids ~= 0);
    jdx = elt.var_ids(j);

    % There are nlocals local position variables.  Allocate 2*nlocals
    % space so we can store both the position and the velocity variables.

    nlocals = length(elt.var_ids);
    qlocal = zeros(2*nlocals,1);

    qlocal(j)         = q(jdx);         % Get the position variables
    if (has_xdot)
      qlocal(j+nlocals) = q(jdx+net.dof); % Get the velocity variables
    end

    % If we actually have an x-dependent forcing term for this model,
    % plug the local piece of dF/dx into the global dF/dx

    [NonlinearKlocal] = feval(elt.model, 'K', elt.R, elt.parameter, qlocal, t);

    if (~isempty(NonlinearKlocal))
      NonlinearK(jdx,jdx) = NonlinearK(jdx,jdx) + NonlinearKlocal(j,j);
    end

  end
end
