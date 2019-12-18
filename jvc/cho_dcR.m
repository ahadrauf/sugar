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

function [q, converged] = cho_dcR(net, q0, is_sp)

if (nargin == 2) & (~isempty(q0)) 
  q = q0;
else
  q = zeros(net.dof,1);
end

if (nargin < 3)
  is_sp = 1;
end

[K] = assemble_system(net, 'K', is_sp);

%T
T=eye(net.dof);
for i = 1:length(net.elements)
   Ttemp=eye(net.dof);
   elt = net.elements(i);
   [Tlocal] = feval( elt.model,  'T', elt.R, elt.parameter );
   if (~isempty(Tlocal))
      j = find(elt.var_ids ~= 0);
      jdx = elt.var_ids(j);
      Ttemp(jdx,jdx) = Tlocal(j,j);
      T=Ttemp*T;
   end
end

%T1
T1=eye(net.dof);
for i = 1:length(net.elements)
   elt = net.elements(i);
   [Tlocal] = feval( elt.model,  'T1', elt.R, elt.parameter );
   if (~isempty(Tlocal)) & i==3
      j = find(elt.var_ids ~= 0);
      jdx = elt.var_ids(j);
      T1(jdx,jdx) = Tlocal(j,j);
   end
end

%T2
T2(1:6,1:6)=eye(6);
T2(7:12,1:6)=eye(6);
T2(13:18,13:18)=eye(6);

T=T*T1*T2;

%r=rot2local(0,0,pi/4);

G = K*q - assemble_F(net, q, 0);
dG = K - assemble_dFdx(net, q, 0, is_sp);
delta_q = (dG\G);

iter = 0;

while (norm(delta_q ./ net.scales') >= 1e-6 & iter < 40) 
  q = q - delta_q;
  G = K*q - assemble_F(net, q, 0);
  dG = K - assemble_dFdx(net, q, 0, is_sp);
  delta_q = (dG\G);
  iter = iter + 1;
end

q = q - delta_q;

converged = 1;
if iter == 40
  disp('Warning: DC solution finder did not converge after 40 iterations');
  converged = 0;
end
