%net=cho_load('can6.m');q=cho_dc(net);figure(1);cho_display(net,q);  
 
uses mumps2.net

anchor p1 [A]   [l=10u w=10u oz=-90]
beam3d p1 [A B] [l=100u w=2u h=2u oz=90 ox=0 oy=0]
%freehinge p1 [B C] [l=100u w=20u h=2u oz=45 ox=0 oy=0]
beam3d p1 [B C] [l=100u w=20u h=2u oz=45]
beam3d p1 [C D] [l=100u w=2u h=2u oz=0]
anchor p1 [D]   [l=10u w=10u]

f3d * [C] [F=30u oz=-90]
%f3d * [C] [F=30u oy=90]

