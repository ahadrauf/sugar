% net=cho_load('rigid1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('rigid1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('rigid1.m');figure(1);cho_display(net);

% net=cho_load('rigid1.m');[q]=cho_dcR(net);figure(1);cho_display(net,q);

uses mumpsx.net
anchor p1 [A]   [l=10u w=10u oz=180]
beam3d p1 [A B] [l=200u w=8u h=2u]
beam_R2 p1 [B C] [l=200u w=2u h=2u oz=0]
beam3d p1 [C D] [l=200u w=8u h=2u]
f3d * [D][F=1000000e-6 oz=0] 

