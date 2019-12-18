%net=cho_load('can8.m');q=cho_dc(net);figure(1);cho_display(net,q);  

uses mumps2.net

anchor     p1 [A]   [l=10u w=10u]
hi beam3d  p1 [A B] [l=100u w=6u h=2u oz=0]
bye beam3d p1 [B C] [l=100u w=6u h=2u oz=45 j(1)=5]
f3d * [C] [F=100u oz=-90]

