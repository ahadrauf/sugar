% [q,K,Knl,F,dFdx] = cho_dc_cubic(net, q0, is_sparsearse)    
% Sugar function. Static analysis of a cubic nonlinear stiffness.
% Performs a simple Newton-Raphson iteration to find an equilibrium 
% point of K*q + Knl*q^3= F(q).
% Iteration stops if step size < 1e-9, or after 40 steps.
% 
% Inputs:
%   net - netlist structure
%   q0  - (Optional) - Starting guess of equilibrium position
%         If q0 is not provided, the search will start at q0 = 0.
%   is_sparse - (Optional) - Flags whether to use sparse solvers or not.
%         Default is true (use sparse solvers).
% 
% Outputs:
%   q    - estimated equilibrium state
%   K    - stiffness constant
%   Knl  - cubic nonlinear stiffness constant 
%   F    - last force
%   dFdx - last jacobian 
%
% See also
%   cho_dc

function [q,K,Knl,F,dFdx] = cho_dc_cubic(net, q0, is_sparsearse) 
tolerance=1e-9; %Newton-Raphson tolerance limit
iterationlimit=40; %limit of iterations
true=1; 

q = zeros(net.dof,1); %allocate q0
if (nargin==2)&(~isempty(q0)) q=q0; end %user-defined starting guess
if (nargin<3) is_sparse=true; end %use sparse matrix methods

%initialize the Newton-Raphson loop
[K] = assemble_system(net, 'K', is_sparse); %linear stiffness constant
[Knl]=assemble_system(net, 'Knl', is_sparse); %cubic nonlinear stiffness constant
F = assemble_F(net,q,0); %excitation
dFdx = assemble_dFdx(net,q,0,is_sparse); %jacobian
G = K*q + Knl*q.^3 - F; %value of G(q)
dG = K + 3*Knl*diag(q.^2) - dFdx; %slope of G(q)
delta_q = dG\G; %change in displacement
iteration=0; %reset interation counter

%do Newton-Raphson interations while the change in displacement is within tolerance and iterations < limit.
while ((norm(delta_q)>=tolerance)&(iteration<iterationlimit)) 
   q = q - delta_q; %update displacement
   F = assemble_F(net,q,0); %update force
   dFdx = assemble_dFdx(net,q,0,is_sparse); %update jacobian
   G = K*q + Knl*q.^3 - F; %value of G(q) 
   dG = K + 3*Knl*diag(q.^2) - dFdx; %slope of G(q)
   delta_q = dG\G; %change in displacement
   iteration=iteration+1; %interation count
   
   g=Knl*q.^3 - F;
   g(2)
   
end

q = q - delta_q; %update displacement one last time
if ((iteration==iterationlimit)&(norm(delta_q)>tolerance)) %nonconvergence warning
   fprintf('\nWarning: Equilibrium solution finder did not converge after %d iterations.\n',iterationlimit);
end

