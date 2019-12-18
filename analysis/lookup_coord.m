% Look up the index of a coordinate in the global state vector.
%
% Inputs:
%   net   - netlist structure
%   node  - name of node
%   coord - name of coordinate at node
%
% Output:
%   id    - index of the indicated coordinate.  If the coordinate
%           could not be found, or if it is grounded, then id == 0.

function [id] = lookup_coord(net, node, coord)

for i=1:length(net.nodes)
   if (strcmp(net.nodes(i).name, node))   
      for j=1:length(net.vars)
         if (net.vars(j).owner == i & ~strcmp(net.vars(j).type, 'b') & strcmp(net.vars(j).name, coord))            
            if (j > net.dof)
               id = 0;
            else
               id = j;
            end
            return
         end
     end
     id = [];
     return
   end
end
id = [];

