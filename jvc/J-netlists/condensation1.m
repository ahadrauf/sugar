% net=cho_load('condensation1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('condensation1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('condensation1.m');figure(1);cho_display(net);

uses mumpsx.net
anchor p1 [A]   [l=10u w=10u oz=180]
beam3d p1 [A B] [l=100u w=2u h=2u oz=45]
beam3d p1 [B C] [l=100u w=2u h=2u oz=-45]
beam3d p1 [C D] [l=100u w=2u h=2u oz=0]
f3d * [D][F=2e-6 oz=90] 

