% net=cho_load('can00.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('can00.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('can00.m');figure(1);cho_display(net);

uses mumps.net
anchor p1 [A]   [l=10u w=10u oz=pi]
beam3d p1 [A a] [l=200u w=20u h=2u ]
beam3d p1 [a B] [l=200u w=2u h=2u oz=pi/10]
f3d * [B][F=6e-6 oz=pi] 
