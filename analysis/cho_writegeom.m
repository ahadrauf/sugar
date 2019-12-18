% Writes out a description of the structure, to be used by
% the Java display applet for the web service.
%
% Inputs:
%   net - the netlist structure
%   q   - The displacement of the original structure
%         If q is empty, the undisplaced structure will
%         be written.
%   fname - Name of the file to write

function cho_writegeom(net,q,fname);


% --- First, open the file for write (IEEE fp, big-endian ordering)

fid = fopen(fname, 'w', 'ieee-be');


% --- fid == -1 means just count number of beams

beamCount = 0;
for i = 1 : length(net.elements)

  elt = net.elements(i);

  % TODO: Make intelligent use of varargin.  Also, don't pass in elt.foo,
  %  just pass in elt!
  count = feval( elt.model, 'writebin', elt.R, elt.parameter, [], 0, ...
                 net.nodes(elt.node_ids), -1, []);

  if (~isempty(count))
    beamCount = beamCount + count;
  end

end


% --- write beam info to fid

fwrite(fid, beamCount, 'integer*4');

for i = 1 : length(net.elements)

  elt = net.elements(i);
  feval( elt.model, 'writebin', elt.R, elt.parameter, [], 0, ...
         net.nodes(elt.node_ids), fid, elt.var_ids);

end


% --- Finally, write out the geometry info

fwrite(fid, net.dof, 'integer*4');

if (isempty(q))
  fwrite(fid, zeros(net.dof,1), 'real*4');
else
  fwrite(fid, q, 'real*4');
end

fclose(fid);

