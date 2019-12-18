%net=cho_load('can4.m');q=cho_dc(net);figure(1);cho_display(net,q);  
%f=0.0005e-6;i=1;k(i,1)=f/q(lookup_coord(net,'F','x'));k(i,2)=f/q(lookup_coord(net,'F','y'));k(i,3)=f/q(lookup_coord(net,'F','rz')) 
uses mumps2.net

anchor p1 [A]   [l=10u w=10u]
beam3d p1 [A B] [l=100u w=2u h=2u oz=0]
beam3d p1 [B C] [l=100u w=2u h=2u oz=-90]
beam3d p1 [C D] [l=100u w=2u h=2u oz=-45]
beam3d p1 [D E] [l=150u w=2u h=2u oz=90]
beam3d p1 [E F] [l=50u w=2u h=2u oz=45]


f3d * [F] [F=5u oz=0]

