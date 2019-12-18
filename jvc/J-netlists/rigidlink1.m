% net=cho_load('rigidlink1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'C','y'))
% net=cho_load('rigidlink1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('rigidlink1.m');figure(1);cho_display(net);

uses mumpsx.net
a=pi/2
anchor p1 [A]   [l=10u w=10u oz=-a]
beam3d p1 [A B] [l=100u w=20u h=2u oz=a]
beam3d p1 [B C] [l=100u w=2u L1=50u L2=10u oz=0]
%rigidlinkbeam p1 [B C] [l=90u w=2u L1=10u L2=0u oz=0]
f3d * [C][F=5e-6 oz=a] 
