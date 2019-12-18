%net=cho_load('joint2.m');q=cho_dc(net);figure(1);cho_display(net,q);

uses mumps2.net

anchor  p1 [A]   [l=10u w=10u oz=-(90+45)]
beam3d  p1 [A B] [l=100u w=2u h=2u oz=90]
beam3d  p1 [B C] [l=50u w=2u h=2u oz=-(90+35)]
beam3d  p1 [C D] [l=50u w=2u h=2u oz=-(90+35)]
hinge   p1 [D]   [l=10u w=10u oz=-(90+35)]
%freejoint   p1 [B]   [l=10u w=10u]
%f3d     *  [C]   [F=700u oz=-90] 

