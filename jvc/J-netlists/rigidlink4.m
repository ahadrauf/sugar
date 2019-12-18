% net=cho_load('rigidlink4.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'C','rz')),q(lookup_coord(net,'D','rz'))
% net=cho_load('rigidlink4.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('rigidlink4.m');figure(1);cho_display(net);

uses mumpsx.net
a=pi/2
anchor p1 [A]   [l=10u w=10u oz=-a]
beam3d p1 [A B] [l=40u w=20u h=2u oz=pi/2]

rigidlinkbeamcorner p1 [B C] [l=100u w=10u L1=10u oz1=0 L2=10u oz2=0 oz=pi/4*0]
%rigidlinkbeamcloaked p1 [B C] [l=100u w=10u L1=100u oz1=pi/4*0 L2=10u oz2=0 oy=pi/2]

beam3d p1 [C D] [l=20u w=2u oz=pi/4]

f3d * [C][F=-50e-6 oz=pi/2] 
%f3d * [D][F=20e-6 oz=0] 

