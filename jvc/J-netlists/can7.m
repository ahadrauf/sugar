%net=cho_load('can7.m');q=cho_dc(net);figure(1);cho_display(net,q);  
 
uses mumps2.net

anchor p1 [A]   [l=10u w=10u]
beam3d p1 [A B] [l=100u w=2u h=2u oz=90]
freehinge p1 [B C] [l=100u w=40u h=2u oz=45 ox=0 oy=0]
beam3d p1 [C D] [l=100u w=2u h=2u oz=0]
beam3d p1 [D E] [l=100u w=2u h=2u oz=0]
anchor p1 [E]   [l=10u w=10u]

freehinge p1 [D F] [l=100u w=40u h=2u oz=45 ox=0 oy=0]
beam3d p1 [F G] [l=100u w=2u h=2u oz=0]
anchor p1 [G]   [l=10u w=10u]

f3d * [D] [F=30u oz=90]

