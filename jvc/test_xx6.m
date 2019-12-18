%slider

K=...
   [1  0  0     0
    0  1  0    -1
    0  0  0.64 -0.64
    0 -1 -0.64  1.64]
 
 s=sqrt(2)/2;
 
T=...
   [s s 0 0
   -s s 0 0
    0 0 1 0
    0 0 0 1]
 
A=...
  [1 0 0 0
   0 0 1 0
   0 0 0 1]'

%example 6.5-Sack


%assemble K
%assemble T
%assemble A
f=4e4; 
Kp=T'*f*K*T
Kpa=A'*Kp*A
%Newton-Raphson loop
%assemble P
P=[0 0 90 0]'
Pp=T'*P
Ppa=A'*Pp
Upa=Kpa\Ppa
U=T*A*Upa

f=4e4; 
TA=T*A
Kpa=TA'*f*K*TA
P=[0 0 90 0]'
Ppa=TA'*P
Upa=Kpa\Ppa
U=TA*Upa

