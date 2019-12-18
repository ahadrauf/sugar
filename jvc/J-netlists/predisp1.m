%net=cho_load('predisp1.m');figure(1);cho_display(net);
%net=cho_load('predisp1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net
%uses predisplacedbeam2.m
anchor p1 [a1][l=10u h=10u w=20u oz=pi]
predisplacedbeam4 p1 [a1 a11][l=100u w=10u h=2u qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=55u qz2=0 qox2=0 qoy2=0 qoz2=0] 

f3d * [a11] [F=300u oy=pi/2]