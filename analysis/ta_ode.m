% Evaluates the right hand side of the ODE
%  C [qhat; qhat']' = -G [qhat; qhat'] + D1*F(q)
% where D1 is a diagonal scaling matrix and qhat = D2*q is
% a diagonal scaling of the state vector.
%
% When the flag is 'mass', returns the "mass" matrix C.  Note
% that this is not the same as the original mass matrix M
% in the second order system
%  M q'' + D q' + K q = F(q)
%
% When the flag is 'jacobian', returns the Jacobian of the right
% hand side of the ODE (with respect to [qhat; qhat']):
%  -G + [D1*(dF/dq)*Dq 0; 0 0]
%
% Used internally by the transient simulation routine cho_ta.

function [qdot] = ta_ode(t,q,flag,net,C,G,d1,d2)

if strcmp(flag, 'mass')

  qdot = C; %Store MASS into solver, so we only return RHS.

elseif strcmp(flag, 'jacobian')
   
   DF1 = assemble_dFdx(net, d2.*q, t);
   DF2 = assemble_dFdxdot(net, d2.*q, t);
   J=-G;
   J(1:net.dof, 1:net.dof) = J(1:net.dof, 1:net.dof) + diag(d1(1:net.dof)) * DF1 * diag(d2(1:net.dof)); 
   J(1:net.dof, net.dof+1:end) = J(1:net.dof, net.dof+1:end) + diag(d1(1:net.dof)) * DF2 * diag(d2(1:net.dof)); 
   qdot = J;

else
   
  F = assemble_F(net, d2.*q, t);
  D1F = d1(1:length(F)).*F;  % Giving errors - changed by Ahad on
%   12/18/19
%   D1F = d1(1:length(F), 1:length(F))*F;
  
  qdot = -G*q;
  qdot(1:length(D1F)) = qdot(1:length(D1F)) + D1F(1:length(D1F))';
  
end

