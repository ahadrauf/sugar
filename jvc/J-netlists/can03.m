%beam and slider
% net=cho_load('can03.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('can03.m');figure(1);cho_display(net);

uses mumps.net
a=60
b=a*3.1415/180
anchor p1 [C]   [l=10u w=10u oz=-90]
beam3d p1 [A C] [l=200u w=2u h=2u oz=-a]
slider2 p1 [A B] [l=200u w=2u h=2u oz=a slideroz=b]
beam3d p1 [B BB] [l=10u w=10u h=10u oz=a]

f3d * [A][F=-30u] 

beam3d p1 [A1 C]  [l=200u w=1e-22 h=1u oz=-a]
beam3d p1 [A1 B1] [l=50u w=1e-22 h=1u oz=a]
beam3d p1 [B1 B2] [l=200u w=4u h=1u oz=a]
