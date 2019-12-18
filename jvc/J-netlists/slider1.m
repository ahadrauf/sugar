%beam and slider
% net=cho_load('slider1.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('slider1.m');figure(1);cho_display(net);
%for f=0:-5e-6:-40e-6, p.f=f;net=cho_load('slider1.m',p);[q]=cho_dc_slider(net);figure(1);cho_display(net,q); end

uses mumps.net
param f=-30u
a=pi/3
anchor p1 [C]   [l=10u w=10u oz=-pi/2] %lower anchor
beam3d p1 [A C] [l=200u w=2u h=2u oz=-a] %lower beam
slider2 p1 [A B] [l=200u w=2u h=2u oz=a slideroz=a] %slider beam
f3d * [A][F=f] %horizontal force

%f3d * [A][F=-30u] 
%beam3d p1 [C CC] [l=400u w=8u h=2u oy=pi/4] %lower beam
%slider2 p1 [A D] [l=200u w=2u h=2u oz=0 slideroz=-pi/2] %slider beam, slider on 2nd node
%beam3d p1 [B BB] [l=10u w=10u h=10u oz=a] %stub n/a
%beam3d p1 [A1 C]  [l=200u w=1e-22 h=1u oz=-a]
%beam3d p1 [A1 B1] [l=50u w=1e-22 h=1u oz=a]
%beam3d p1 [B1 B2] [l=200u w=4u h=1u oz=a]

