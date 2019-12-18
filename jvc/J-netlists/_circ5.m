%beam and slider
% net=cho_load('_circ5.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('_circ5.m');figure(1);cho_display(net);
% net=cho_load('_circ5.m');cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('_circ5.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);
% net=cho_load('_circ5.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);q(lookup_coord(net,'g','ry'))

uses mumps.net
W=20u
H=2u
anchor   p1 [a]  [l=10u w=10u h=20u oz=-pi/2]
circbeam2 p1 [a b][l=180u w=W h=H radius=100u alpha=pi/3]
circbeam2 p1 [b c][l=180u w=W h=H radius=100u alpha=pi/3 oz=pi/3]
circbeam2 p1 [c g][l=180u w=W h=H radius=100u alpha=pi/3 oz=2*pi/3]
circbeam2 p1 [g e][l=180u w=W h=H radius=100u alpha=pi/3 oz=pi]
circbeam2 p1 [e f][l=180u w=W h=H radius=100u alpha=pi/3 oz=4*pi/3]
circbeam2 p1 [f d][l=180u w=W h=H radius=100u alpha=pi/3 oz=5*pi/3]
f3d * [d][F=0*500u oy=-pi/2]


