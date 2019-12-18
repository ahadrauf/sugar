% net=cho_load('beamplate2.m');q=cho_dc(net);figure(1);cho_display(net,q);view(0,0)
%lookup_coord(net,'b','z')

uses mumps2.net

anchor p1 [A][l=10u w=10u h=10u oz=180] 
anchor p1 [B][l=10u w=10u h=10u oz=180] 

beam3d p1 [A a][l=200u w=100u h=2u oz=0]
beam3d p1 [B b][l=200u w=100u h=2u oz=0]
beam3d p1 [b a][l=200u w=100u h=2u oz=90]

f3d * [a][ F=400u oy=-90] 
f3d * [b][ F=400u oy=-90] 

