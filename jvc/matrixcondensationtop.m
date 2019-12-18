%This transforms k*q to khat*q1. 
function [khat]=matrixcondensationtop(k)
s=length(k);
%partition matrix
k11=k(1:12,1:12); 
k12=k(1:12,13:s); 
k21=k(13:s,1:12);
k22=k(13:s,13:s);
%solving for khat st khat*q1
khat=k11-k12*(k22\k21);

