% net=cho_load('rigid2.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('rigid2.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q
% net=cho_load('rigid2.m');figure(1);cho_display(net);

uses mumpsx.net
anchor p1 [A]   [l=10u w=10u oz=180]
beam3d p1 [A B] [l=200u w=2u h=2u]
beam3d p1 [B C] [l=200u w=2u h=2u]
%beam3d p1 [C D] [l=200u w=2u h=2u]
%f3d * [D][F=1e-6 oz=-90] 

%f3d * [B][F=-1e-6 oz=-90] 
f3d * [C][F=3e-6 oz=-90] 
%f3d * [B][M=-200e-6*1e-6 oy=-90] 

%f3d * [B][F=1e-6 oz=-90] 
%f3d * [B][M=400e-6*1e-6 oy=-90] 
