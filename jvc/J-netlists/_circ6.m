%beam and slider
% net=cho_load('_circ6.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('_circ6.m');figure(1);cho_display(net);
% net=cho_load('_circ6.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('_circ6.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);
% net=cho_load('_circ6.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);q(lookup_coord(net,'g','ry'))
% q=zeros(12,1);q(4)=pi/4;figure(1);cho_display(net,q)

uses mumps.net
W=10u
H=2u
anchor   p1 [d]  [l=10u w=10u h=20u oz=-pi/2]
circbeam6 p1 [a d][w=W h=H radius=100u alpha=pi/4 oz=pi/2]
%circbeam3 p1 [d dd][w=W h=H radius=100u alpha=pi/4 oz=pi/2]

f3d * [a][M=0.910000u oy=-pi/2]
