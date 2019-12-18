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

function [M] = assemble_removal(net, mflag, is_sp);
if nargin<3, is_sp=0; end
if (is_sp), M = sparse(eye(net.dof)); else M = eye(net.dof); end

for i = 1:length(net.elements)
  elt = net.elements(i);
  [v] = feval( elt.model,  mflag, elt.R, elt.parameter );
  if (~isempty(v)) %Only if a model returns an 'v' vector
     j=find(v~=0); %r should be a 1D vector containing nonzeros signifying rows to remove.
     for k=length(j):-1:1 %for each j found
        jdx=elt.var_ids(j(k)); %select row from system size matrix
        M(:,jdx)=[]; %remove rows bottom up.
     end
  end
end

