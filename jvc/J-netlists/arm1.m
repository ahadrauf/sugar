%net=cho_load('arm1.m');q=cho_dc(net);figure(1);cho_display(net,q);  
 
uses mumps2.net

hinge p1 [A]   [l=10u w=10u h=10u]
beam3d p1 [A B] [l=100u w=20u h=20u oz=45]
beam3d p1 [B C] [l=100u w=20u h=20u oz=45]
beam3d p1 [C D] [l=300u w=20u h=20u oz=-55]

freehinge p1 [b B] [l=200u w=10u h=2u oz=20]
beam3d p1 [a1 b] [l=100u w=2u h=2u oy=90]
beam3d p1 [a2 b] [l=100u w=2u h=2u oy=-90]
anchor p1 [a1]   [l=10u w=10u]
anchor p1 [a2]   [l=10u w=10u]

%beam3d p1 [a1 aa1] [l=100u w=200u h=20u oy=90]

%freehinge p1 [D F] [l=100u w=40u h=2u oz=45]
%beam3d p1 [F G] [l=100u w=2u h=2u oz=0]
%anchor p1 [G]   [l=10u w=10u]

f3d * [b] [F=-100u oz=0]

