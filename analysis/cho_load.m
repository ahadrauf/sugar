% Loads a netlist description into SUGAR's internal format.
%
% Inputs: 
%   netname - netlist file name
%   param   - (Optional) parameter settings
%
% Outputs:
%   net - netlist structure
%
% Comments:
%   This function supercedes the old cho_load function, which is
% now available as cho_load_old.  Users are encouraged to port over
% any netlists in the old format to the new format.

function [net] = cho_load(netname, param);

if ~exist('sugar_c')
  disp('The sugar_c MEX file is not in your path, and may not be');
  disp('compiled for your system.  To compile sugar_c for your');
  disp('system, change to the SUGAR compile subdirectory (from');
  disp('within Matlab) and type "makemex".  Then copy the resulting');
  disp('mex file (sugar_c.mex??? or sugar_c.dll) to the analysis');
  disp('subdirectory.'); 
  error('sugar_c does not exist');
end

if (nargin == 1)
  net = sugar_c(netname);
else
  net = sugar_c(netname, param);
end


% -- Position post-processing

for i = 1:length(net.elements)

elt = net.elements(i);

  output = feval( elt.model, 'postpos', elt.R, elt.parameter, [], [], ...
                 net.nodes(elt.node_ids) );

  if (isfield(output, 'params'))
    net.elements(i).parameter = output.params;
  end


  if (isfield(output, 'R'))
    net.elements(i).R = output.R;
  end

end

% -- Scale post-processing

for i = 1:net.dof
  net.scales(i) = default_scales(net.vars(i).name);
end

