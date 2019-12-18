%This transforms [q1;q2;q3] to [q1;q2] where q1 and q2 are the end nodes of the circringpos1
function [khat]=matrixcondensation_2endnodecirc(k)
s=length(k);
%break k into [k1,k2,k3;k4,k5,k6;k7,k8,k9]
k11=k(1:12,1:12);
k12=k(1:12,13:s);
k21=k(13:s,1:12);
k22=k(13:s,13:s);
%condensed matrix into a 12x12, corresponding to [q1;q2].
khat=k11-k12*(k22\k21);


