%beam and slider
% net=cho_load('_circ1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('nettech3.m');figure(1);cho_display(net);

uses mumps.net
param q1=4


anchor   p1 [a]  [l=1u w=1u h=1u]
circbeam p1 [a b][w=10u h=2u radius=100u alpha=pi/3 q=q1(1)]
f3d * [b][F=1000u ]

