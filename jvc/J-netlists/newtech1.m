%beam and slider
% net=cho_load('newtech3.m');figure(1);cho_display(net);

uses mumps.net
param qax=0,qay=0,qaz=0,qaox=0,qaoy=0,qaoz=0
param qbx=0,qby=0,qbz=0,qbox=0,qboy=0,qboz=0
param qcx=0,qcy=0,qcz=0,qcox=0,qcoy=0,qcoz=0
param F=0, Foz=0, M=0

L=300u
W=2u
H=2u

anchor   p1 [A]  [l=10u w=10u h=10u]
%predisplacedbeam4 p1 [A a][l=L w=W h=H qox1=0    qoy1=0    qoz1=0      qx2=qax qy2=qay qz2=qaz qox2=qaox qoy2=qaoy qoz2=qaoz ] 
beam3d p1 [A a][l=L w=W h=H qox1=0    qoy1=0    qoz1=0      qx2=qax qy2=qay qz2=qaz qox2=qaox qoy2=qaoy qoz2=qaoz ] 

f3d * [a][F=F oz=Foz]
f3d * [a][M=M oy=-pi/2]

