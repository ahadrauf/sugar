%net=cho_load('joint3.m');q=cho_dc(net);figure(1);cho_display(net,q);

uses mumps2.net

anchor  p1 [A]   [l=10u w=10u oz=-90]
beam3d  p1 [A B] [l=100u w=2u h=2u oz=90]
beam3d  p1 [B M] [l=100u w=2u h=2u oz=0]
beam3d  p1 [M C] [l=100u w=2u h=2u oz=0]
beam3dh  p1 [C D] [l=100u w=2u h=2u oz=-90]
hinge   p1 [D]   [l=10u w=10u oz=-90]
%freehinge   p1 [C]   [l=10u w=10u]
%f3d     *  [M]   [F=100u oz=90] 

