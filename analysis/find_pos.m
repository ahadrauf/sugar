%find_pos
% Locate the mechanical nodes in the system and determine their position.
%
% Inputs:
%  elements, nodes -- the element and node lists from the netlist
%
% Outputs:
%  nodes augmented with .pos fields indicating undisplaced position

% David Bindel, 5/26/00; 
%  modified 8/26 to work with experimental 2.0 structs;
%  modified 9/1 to move relative position info into model functions
%  modified 11/11 to permit models with mixed mechanical and non-mech
%    nodes (including branch nodes).  Note that mechanical nodes must
%    appear first.
%  modified 11/11 to allow for absolute positions of nodes.

function [nodes] = find_pos(elements, nodes);

% Locate mechanical elements / nodes

j = 0; 
mech_elts = zeros(length(elements));
mech_nodes = [];
q = [];
nodepos = zeros(length(nodes), 3);

for i = 1 : length(elements)
  elt = elements(i);

  relpos = feval(elt.model,  'pos', elt.R, elt.parameter);
  if ~isempty(relpos)
    [ncoords,nnodes] = size(relpos);
    mech_elts(i) = 1;
    mech_nodes = union( mech_nodes, elt.node_ids(1:nnodes) );
  end

  abspos = feval(elt.model,  'abspos', elt.R, elt.parameter);
  if ~isempty(abspos)
    [ncoords,nnodes] = size(abspos);
    mech_elts(i) = 1;
    mech_nodes = union( mech_nodes, elt.node_ids(1:nnodes) );
    q = union( q, elt.node_ids(1:nnodes) );
    nodepos(elt.node_ids(1:nnodes), :) = abspos'; 
  end

end

if (isempty(q) & ~isempty(mech_nodes))
  q = mech_nodes(1);
end

% Do a tree search to compute relative coordinates of nodes.

nodes_left = setdiff(mech_nodes, q);
nodecolor = ones(length(nodes));
nodecolor(nodes_left) = 0;

while (~isempty(nodes_left))

  if (isempty(q))
    disp('Warning: Not all mechanical nodes are connected!');
    current_node = nodes_left(1);
    nodecolor(current_node) = 1;
    q = [ current_node ];
  end

  while (~isempty(q))

    current_node = q(1);
    current_pos = nodepos(current_node,:)';
    q = q(2:end);

    % Go through each mechanical element and align node positions
    for elt_idx = nodes(current_node).elt_ids
      if mech_elts(elt_idx)

        elt = elements(elt_idx);
        uncolored = elt.node_ids(find(nodecolor(elt.node_ids) == 0));

        if ~isempty(uncolored)

          % Color unmarked children and add them to the queue
          q = union(q, uncolored);
          nodecolor(uncolored) = 1;

          % Figure out the relative positions of all the nodes in this elt
          relpos = feval(elt.model, 'pos', elt.R, elt.parameter);
          [ncoords,nnodes] = size(relpos);

          % Translate the relative positions so that the current node
          % is at its current position. 
          current_node_idx = find(elt.node_ids == current_node);
          offsets = repmat(current_pos-relpos(:,current_node_idx), 1, ...
                           nnodes);
          relpos = relpos + offsets;
          nodepos(elt.node_ids(1:nnodes), :) = relpos';

        end
      end 
    end
    
  end

  nodes_left = find(nodecolor == 0);
end

% Add the node positions into the node information structure
for i = mech_nodes
  nodes(i).pos = nodepos(i,:)';
end
