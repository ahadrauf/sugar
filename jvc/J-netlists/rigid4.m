% net=cho_load('rigid4.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('rigid4.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('rigid4.m');figure(1);cho_display(net);

uses mumpsx.net
anchor p1 [A]   [l=20u  w=20u oz=180]
beam3d p1 [A B] [l=100u w=2u h=2u oz=45]
beam3d p1 [B C] [l=100u w=8u h=2u oz=-45]
beam3d p1 [C D] [l=200u w=2u h=2u oz=0]
beam3d p1 [D E] [l=200u w=8u h=2u oz=90]
beam3d p1 [E F] [l=200u w=2u h=2u oz=90]

f3d * [F][F=1e-6 oz=0] 

