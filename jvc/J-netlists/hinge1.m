% net=cho_load('hinge1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('hinge1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('hinge1.m');figure(1);cho_display(net);

uses mumpsx.net
anchor p1 [A]   [l=10u w=10u oz=pi]
beam3d p1 [A B] [l=200u w=2u h=2u]
hinge2 p1 [B C] [l=210u w=8u h=2u oz=0]
beam3d p1 [C D] [l=200u w=2u h=2u oz=0]
beam3d p1 [D E] [l=200u w=2u h=2u oz=-pi/2]
anchor p1 [E]   [l=10u w=10u oz=0]
f3d * [C][F=100000e-6 oz=0] 

