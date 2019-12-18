function [qdot] = ta_ode_regular(t,q,flag,net,M,D,K,l,u,p,il,iu,ip,invMD,invMK)

%if strcmp(flag, 'mass')
%
%  qdot = C; %Store MASS into solver, so we only return RHS.
%
%elseif strcmp(flag, 'jacobian')
%   
%  qdot = G;
%  qdot(1:net.dof, 1:net.dof) = qdot(1:net.dof, 1:net.dof) + ... 
%                               assemble_dF(net, q, t) * ...
%                                 diag(d2(1:net.dof));
%
%else
   
   
  F = assemble_F(net, q, t);
  
%  qdot = [q(net.dof+1:2*net.dof); M\(F - D*q(net.dof+1:2*net.dof) - K*q(1:net.dof))];
%  qdot = [q(net.dof+1:2*net.dof); (M\F - invMD*q(net.dof+1:2*net.dof) - invMK*q(1:net.dof))];
%  qdot = [q(net.dof+1:2*net.dof); ( u\(l\(p*F)) - invMD*q(net.dof+1:2*net.dof) - invMK*q(1:net.dof))];
%  qdot = [q(net.dof+1:2*net.dof); iu*(il*(p*(F - D*q(net.dof+1:2*net.dof) - K*q(1:net.dof))))];
%  qdot = [q(net.dof+1:2*net.dof); iu*il*p*(F - D*q(net.dof+1:2*net.dof) - K*q(1:net.dof))];
qdot = [q(net.dof+1:2*net.dof); u\(l\(p*(F - D*q(net.dof+1:2*net.dof) - K*q(1:net.dof))))];


%  qdot = [q(net.dof+1:2*net.dof); (((u\l)\p)*(F - D*q(net.dof+1:2*net.dof) - K*q(1:net.dof)))];
%  qdot = [q(net.dof+1:2*net.dof); u\l\p*(((F - D*q(net.dof+1:2*net.dof) - K*q(1:net.dof))))];
%  qdot = [q(net.dof+1:2*net.dof); p*il*iu*(((F - D*q(net.dof+1:2*net.dof) - K*q(1:net.dof))))];
  
%end

