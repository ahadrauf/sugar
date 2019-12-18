function [q,khat,fhat]=matrixcondensation_endnode(k,force,bottom)
if bottom %if bottom partition of the matrix
   s=length(k);f=sparse(zeros(s,1));
   if length(force)==6
      f(s-6+1:s,1)=force;
   else
      f=force;
   end
   %partition k=[k11,k12;k21,k22] f=[f1;f2] where f2 is the 6dof end node, f1 is the N-6 dof other nodes.
   k11=sparse(k(1:s-6,1:s-6));
   k12=sparse(k(1:s-6,s-6+1:s));
   k21=k12';
   k22=sparse(k(s-6+1:s,s-6+1:s));
   f2=sparse(f(s-6+1:s,1));
   f1=sparse(f(1:s-6,1));
   %condensed k and f
   khat=k22-k21*(k11\k12);
   fhat=f2-k21*(k11\f1);
   %displacement of end node
   q=khat\fhat;
else %if top partition of the matrix
   s=length(k);f=sparse(zeros(s,1));
   if length(force)==6
      f(1:6,1)=force;
   else
      f=force;
   end
   %partition k=[k11,k12;k21,k22] f=[f1;f2] where f1 is the 6dof end node, f2 is the N-6 dof other nodes.
   k11=sparse(k(1:6,1:6));
   k12=sparse(k(1:6,6+1:s));
   k21=k12';
   k22=sparse(k(6+1:s,6+1:s));
   f2=sparse(f(6+1:s,1));
   f1=sparse(f(1:6,1));
   %condensed k and f
   khat=k11-k12*(k22\k21);
   fhat=f1-k12*(k22\f2);
   %displacement of end node
   q=khat\fhat;
end