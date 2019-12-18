function [Z,chi2omega0,omega02,alpha,ZTF,F]=decouplesystem(net)
is_sparse=0; %since not using eigs
K=assemble_system(net,'K',is_sparse); %linear stiffness constant
Knl=assemble_system(net,'Knl',is_sparse); %cubic nonlinear stiffness constant
F=assemble_F(net,zeros(length(K),1),0); %initial excitation
D=assemble_system(net,'D',is_sparse); %linear stiffness constant
M=assemble_system(net,'M',is_sparse); %linear stiffness constant
f.disp=0;
f.tol=1e-120;
%[X,omega02]=eigs(K,M,f);
[X,omega02]=eig(K,M);
N = length(M);
Z=[];
for i = 1 : N
   M(i) = X(1:N,i)'*M*X(1:N,i);
   Zi = (X(1:N,i)/sqrt(M(i)));
   Z = [Z,Zi];
end

omega02=[];
for i = 1 : N
   omega02(i)=Z(1:N,i)'*K*Z(1:N,i); %just a check
   alpha(i)=Z(1:N,i)'*Knl*Z(1:N,i); %just a check
   chi2omega0(i)=Z(1:N,i)'*D*Z(1:N,i);
   I(i)=Z(1:N,i)'*M*Z(1:N,i); %just a check
   ZTF(i)=Z(1:N,i)'*F;
end

F=Z(1:N,2)'*F;
kb=alpha(2)
ka=omega02(2)

ql=F/ka;
ql=Z*ql;ql=ql(2,2)

qnl=1/6*((108*F+12*sqrt(3)*sqrt((4*ka^3+27*F^2*kb)/kb))*kb^2)^(1/3)/kb-2*ka/(((108*F+12*sqrt(3)*sqrt((4*ka^3+27*F^2*kb)/kb))*kb^2)^(1/3));
qnl=Z*qnl;qnl=qnl(2,2)

