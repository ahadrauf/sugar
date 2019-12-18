function [q,K,Tv,F,T,v] = cho_dc_slider(net, q0, is_sp)
if (nargin == 2) & (~isempty(q0)), q = q0; else q = zeros(net.dof,1); end
if (nargin < 3), is_sp = 1; end
is_sp=0;

K = assemble_system(net, 'K', is_sp);
T=assemble_system_from_eye(net,'T',is_sp); %transformation
v=assemble_removal(net,'v'); %removal
Tv=T*v;

%TA=T*A
%Kpa=Tv'*K*Tv
%P=[0 0 90 0]'
%Ppa=Tv'*P
%Upa=Kpa\Ppa
%U=Tv*Upa

G = K*q - assemble_F(net, q, 0);
dG = K - assemble_dFdx(net, q, 0, is_sp);

delta_q_prime = (Tv'*dG*Tv)\(Tv'*G); %transform
delta_q = Tv*delta_q_prime; %convert back

U=(Tv'*K*Tv)\(Tv'*assemble_F(net, q, 0));
U2=K\assemble_F(net, q, 0);

iter = 0;

%Newton-Raphson
while (norm(delta_q) >= 1e-9 & iter < 40)
   q = q - delta_q;
   
   F=assemble_F(net, q, 0);
   G = K*q - F;
   dG = K - assemble_dFdx(net, q, 0, is_sp);
   
   delta_q_prime = (Tv'*dG*Tv)\(Tv'*G); %transform
   delta_q = Tv*delta_q_prime; %convert back
   
   iter = iter + 1;
end


q = q - delta_q;

if iter == 40
  disp('Warning: DC solution finder did not converge after 40 iterations');
end
