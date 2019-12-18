function [q,K,F] = cho_dc_linear(net, q0)
if (nargin == 2) & (~isempty(q0)) 
  q = q0;
else
  q = zeros(net.dof,1);
end
K = assemble_system(net, 'K', 1);
F = assemble_F(net, q, 0); %gravitational force
q = K\F;
