%Simple example: 
[T,Q]=ode45('SecondOderODE',[0;1e-6],[0;0],[],M,D,K) 
%where the ode might be
function qdot=SecondOrderODE(t,q,flag,M,D,K)
%M mass, D damping, K stiffness, F excitation
qdot=[q(2); M\(F(t,q)-D*q(2)-K*q(1))];



