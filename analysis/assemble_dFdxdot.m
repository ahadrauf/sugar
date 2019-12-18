% Assemble the Jacobian dF/dx of the forcing function F
%
% Inputs:
%   net   - netlist structure
%   q     - state vector of [x;xdot]
%   t     - time
%   is_sp - (optional) flag whether the Jacobian should be in sparse form.
%           default is false.
%
% Output:
%   dFdxdot  - Jacobian matrix with respect to x'

function [dFdxdot] = assemble_dFdxdot(net, q, t, is_sp);

if (nargin < 4)
  is_sp = 0;
end

if (is_sp)
  % Go to sparse form... pure guess work on nnz
  dFdxdot = spalloc(net.dof, net.dof, 6*net.dof);
else
  dFdxdot = zeros(net.dof);
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

    [dFdxdotlocal] = feval(elt.model, 'dFdxdot', elt.R, ...
                           elt.parameter, qlocal, t);

    if (~isempty(dFdxdotlocal))
      dFdxdot(jdx,jdx) = dFdxdot(jdx,jdx) + dFdxdotlocal(j,j);
    end

  end
end
