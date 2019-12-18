%beam and slider
% net=cho_load('can05.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('can05.m');figure(1);cho_display(net);

uses mumps.net

pi=3.1416
anchor   p1 [a]  [l=1u w=1u h=1u oz=-90]
circbeam p1 [a b][l=180u w=100u h=10u radius=100u alpha=pi/3]
circbeam p1 [b c][l=180u w=100u h=10u radius=100u alpha=pi/3 oz=60]
circbeam p1 [c d][l=180u w=100u h=10u radius=100u alpha=pi/3 oz=120]
