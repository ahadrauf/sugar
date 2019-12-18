% Do a simple Newton-Raphson iteration to find an equilibrium point
% Kq = F(q).  The iteration will stop after forty steps, or after the
% step sizes decrease to a size less than 1e-9.
%
% Inputs: 
%   net - netlist structure
%   q0  - (Optional) - Starting guess of equilibrium position
%         If q0 is not provided, the search will start at q0 = 0.
%   is_sp - (Optional) - Flags whether to use sparse solvers or not.
%         Default is true (use sparse solvers).
%
% Outputs:
%   q         - the estimated equilibrium state
%   converged - true if the convergence test passed

function [q,K,F] = cho_dckf(net, q0, is_sp)

if (nargin == 2) & (~isempty(q0)) 
  q = q0;
else
  q = zeros(net.dof,1);
end

if (nargin < 3)
  is_sp = 1;
end

[K] = assemble_system(net, 'K', 0);
F=assemble_F(net, q, 0);

G = K*q - assemble_F(net, q, 0);
dG = K - assemble_dFdx(net, q, 0, is_sp);
delta_q = dG\G;

iter = 0;

while (norm(delta_q ./ net.scales') >= 1e-6 & iter < 40)

  q = q - delta_q;
  G = K*q - assemble_F(net, q, 0);
  dG = K - assemble_dFdx(net, q, 0, is_sp);
  delta_q = dG\G;
  iter = iter + 1;

end

q = q - delta_q;

converged = 1;
if iter == 40
  disp('Warning: DC solution finder did not converge after 40 iterations');
  converged = 0;
end
