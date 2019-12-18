%This transforms [q1;q2;q3] to [q1;q3;q2] where q1 and q3 are the end nodes.
function [khat]=matrixcondensation_2endnode(k)
s=length(k);
%break k into [k1,k2,k3;k4,k5,k6;k7,k8,k9]
k1=k(1:6,1:6);
k2=k(1:6,7:s-6);
k3=k(1:6,s-5:s);
k4=k(7:s-6,1:6);
k5=k(7:s-6,7:s-6);
k6=k(7:s-6,s-5:s);
k7=k(s-5:s,1:6);
k8=k(s-5:s,7:s-6);
k9=k(s-5:s,s-5:s);
%transform matrix so that it corresponds to [q1;q3;q2]. And partition matrix into 4 sections [k11,k12;k21,k22].
k11=[k1,k3;k7,k9];
k12=[k2;k8];
k21=[k4,k6];
k22=[k5];
%condensed matrix into a 12x12, corresponding to [q1;q3].
khat=k11-k12*(k22\k21);


