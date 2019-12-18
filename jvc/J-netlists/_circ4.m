%beam and slider
% net=cho_load('_circ4.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('_circ4.m');figure(1);cho_display(net);
% net=cho_load('_circ4.m');cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('_circ4.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);
% net=cho_load('_circ4.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);q(lookup_coord(net,'g','ry'))
% q=zeros(12,1);q(4)=pi/4;figure(1);cho_display(net,q)

uses mumps.net
W=20u
H=2u
anchor   p1 [a]  [l=10u w=10u h=20u oz=-pi/2]
circbeam2 p1 [a d][w=W h=H radius=100u alpha=pi/2 oz=0]
%circbeam2 p1 [d dd][w=W h=H radius=100u alpha=pi/4 oz=pi/4]
f3d * [d][F=0*10u oy=-pi/2]
