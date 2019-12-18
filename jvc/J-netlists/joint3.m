% net=cho_load('joint3.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'C','y'))
% net=cho_load('joint3.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('joint3.m');figure(1);cho_display(net);

uses mumpsx.net
anchor p1 [A]   [l=10u w=10u oz=pi] %left anchor
beam3d p1 [A B] [l=200u w=4u h=2u oz=0] %horizontal beam

%beamjoint12 p1 [B C] [l=200u w=4u h=2u oz=pi/4] %diagonal jointed beam
beam3d p1 [B C] [l=200u w=4u h=2u oz=pi/4] %diagonal jointed beam

beam3d p1 [C D] [l=200u w=4u h=2u oz=-pi/2] %vertical beam
anchor p1 [D]   [l=10u w=10u oz=0] %right anchor
f3d * [B][F=50e-6 oz=pi/2] %vertical force

