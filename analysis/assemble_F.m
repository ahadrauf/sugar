% Assemble the the forcing function F
%
% Inputs:
%   net - netlist structure
%   q   - state vector of [x;xdot] 
%         (note: xdot can be ommitted, in which case it is presumed 0)
%   t   - time
%
% Output:
%   F   - force value 

function [F] = assemble_F(net, q, t);
%Modification Aug2003 - jvclark

F = zeros(net.dof,1);
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

    [Flocal] = feval(elt.model, 'F', elt.R, elt.parameter, qlocal, t, [], i);

    if (~isempty(Flocal))
      F(jdx) = F(jdx) + Flocal(j);
    end
    
 end
end

