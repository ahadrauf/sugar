%net=cho_load('t1');q=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net

anchor p1 [a][l=10u w=10u h=10u]

beam3d p1 [a b][l=100u w=2u oz=0]
%beam3d p1 [a b1][l=50u w=10u oz=pi/4]
%beam3d p1 [a b2][l=50u w=10u oz=pi/4]
%beam3d p1 [a b3][l=50u w=10u oz=pi/4]

f3d * [b][M=0.01u oy=-pi/2]
%f3d * [b1][F=50000u oz=(pi/2+pi/4)]
%f3d * [b2][F=25000u oz=(pi/2+pi/4)]

