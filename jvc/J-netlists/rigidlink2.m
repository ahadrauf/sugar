% net=cho_load('rigidlink2.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'C','rz')),q(lookup_coord(net,'D','rz'))
% net=cho_load('rigidlink2.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('rigidlink2.m');figure(1);cho_display(net);

uses mumpsx.net
a=pi/2
anchor p1 [A]   [l=10u w=10u oz=-a]
beam3d p1 [B C] [l=100u w=2u h=2u oz=0]
beam3d p1 [A B] [l=50u w=2u h=2u oz=0]
anchor p1 [B]   [l=10u w=10u oz=-a]
rigidlinkbeam p1 [A D] [l=100u w=2u L1=0u L2=50u oz2=-pi/2 oz=0]
%beam3d p1 [D A] [l=100u w=2u L1=0u L2=50u oz=pi]
f3d * [C][F=20e-6 oz=0] 
f3d * [D][F=20e-6 oz=0] 

