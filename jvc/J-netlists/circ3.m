%net=cho_load('circ3.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%q(lookup_coord(net,'b','z'));q50=q(1:6,1)
   
uses mumps.net

anchor p1 [a][l=10u w=10u h=10u oz=pi]
semicircularbeam p1 [a b] [radius=100u w=20u h=2u alpha=pi/4 L1=20u]
semicircularbeam p1 [b c] [radius=100u w=20u h=2u alpha=pi/4 oz=pi/4]
anchor p1 [c][l=10u w=10u h=10u oz=pi/2]

%semicircularbeam p1 [a e] [radius=100u w=20u h=2u alpha=pi/3]

f3d * [b][F=3000u oy=pi/2]
%f3d * [e][F=  -500u oy=pi/2]


