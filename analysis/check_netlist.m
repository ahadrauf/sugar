% Check netlist for parameter errors (parameters missing/out of range).
% If errors are detected, diagnostic messages are printed, and the routine
% will return errflag = 1.
%
% Inputs:
%   net - netlist structure to be checked
%
% Outputs:
%   errflag - logical flag indicating whether an error was encountered

function [errflag] = check_netlist(net);

errflag = 0;
for i = 1:length(net.elements)

  elt = net.elements(i);
  if ( ~exist(elt.model, 'file') )
    err = 'Model cannot be found.';
  else
    err = feval( elt.model,  'check', elt.R, elt.parameter);
  end

  if (~isempty(err))
    name = elt.name;
    if (isempty(name))
      name = '[no name]';
    end
    disp(sprintf('In element %s (%s): %s', ...
                 name, elt.model, err));
    errflag = 1;
  end

end
