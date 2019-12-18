%net=cho_load('joint4.m');q=cho_dc(net);figure(1);cho_display(net,q);

uses mumps2.net

anchor  p1 [A]   [l=10u w=10u oz=-90]
beam3d  p1 [A B] [l=100u w=6u h=2u oz=0]
beam3d  p1 [B M] [l=100u w=6u h=2u oz=0]
beam3d  p1 [M C] [l=100u w=6u h=2u oz=0]
hinge   p1 [M]   [l=10u w=5u oz=-90]
f3d     *  [B]   [F=1900u oz=90] 

