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

function [q,K,f, converged] = cho_dc_secant(net, q0, is_sp)

maxiterations=80;

%Get the first guess at the solution
if (nargin == 2) & (~isempty(q0)) 
  q1 = q0;
else
  q1 = zeros(net.dof,1);
end

if (nargin < 3)
  is_sp = 1;
end

%Use a simple Newton method to determine the second guess
f=assemble_F(net, q1, 0);
[K] = assemble_NonlinearK(net, q1, 0, is_sp);
G = K*q1 - assemble_F(net, q1, 0);
dG = K - assemble_dFdx(net, q1, 0, is_sp);
delta_q = dG\G;
q2 = q1 - delta_q/1e3;

%Use the secant method (~order 1.6, no Jacobian)
iter = 0;
while (norm(delta_q ./ net.scales') >= 1e-6 & iter < maxiterations)
    [K] = assemble_NonlinearK(net, q1, 0, is_sp);
    G1 = K*q1 - assemble_F(net, q1, 0);
    [K] = assemble_NonlinearK(net, q2, 0, is_sp);
    G2 = K*q2 - assemble_F(net, q2, 0);
    delta_q = ((G2-G1)\(q2-q1))*G2;
    q3 = q2 - delta_q;
    q1=q2;
    q2=q3;
    iter = iter + 1;
end

q = q3;
iter
converged = 1;
if iter == maxiterations
  disp('__________Warning: DC solution finder did not converge after max iterations._____________');
  converged = 0;
end
