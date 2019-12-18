%net=cho_load('circ5.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%q(lookup_coord(net,'b','z'));q50=q(1:6,1)
   
uses mumps.net
w=198u
anchor p1 [A][l=10u w=10u h=10u oz=pi]
beam3d p1 [A b][l=150u w=20u oz=pi/2]

semicircularbeam p1 [b b1] [radius=100u w=w h=2u alpha=pi/4 oz=0*pi/4]
semicircularbeam p1 [b1 b2] [radius=100u w=w h=2u alpha=pi/4 oz=1*pi/4]
semicircularbeam p1 [b2 b3] [radius=100u w=w h=2u alpha=pi/4 oz=2*pi/4]
semicircularbeam p1 [b3 b4] [radius=100u w=w h=2u alpha=pi/4 oz=3*pi/4]
semicircularbeam p1 [b4 b5] [radius=100u w=w h=2u alpha=pi/4 oz=4*pi/4]
semicircularbeam p1 [b5 b6] [radius=100u w=w h=2u alpha=pi/4 oz=5*pi/4]
semicircularbeam p1 [b6 b7] [radius=100u w=w h=2u alpha=pi/4 oz=6*pi/4]
semicircularbeam p1 [b7 b] [radius=100u w=w h=2u alpha=pi/4 oz=7*pi/4]

beam3d p1 [b4 B][l=150u w=20u oz=pi/2]

f3d * [B][F=10u oy=pi/2]
%f3d * [e][F=  -500u oy=pi/2]


