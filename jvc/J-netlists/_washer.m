%beam and slider
% net=cho_load('_washer.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('_washer.m');figure(1);cho_display(net);

uses mumps.net

pi=3.1416
anchor   p1 [a]  [l=1u w=1u h=1u oz=-90]
circbeam p1 [a b][l=180u w=150u h=10u radius=100u alpha=pi/3]
circbeam p1 [b c][l=180u w=150u h=10u radius=100u alpha=pi/3 oz=60]
circbeam p1 [c d][l=180u w=150u h=10u radius=100u alpha=pi/3 oz=120]
circbeam p1 [d e][l=180u w=150u h=10u radius=100u alpha=pi/3 oz=180]
circbeam p1 [e f][l=180u w=150u h=10u radius=100u alpha=pi/3 oz=240]
circbeam p1 [f g][l=180u w=150u h=10u radius=100u alpha=pi/3 oz=300]

