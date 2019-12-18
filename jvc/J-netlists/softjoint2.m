% net=cho_load('softjoint2.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'C','y'))
% net=cho_load('softjoint2.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('softjoint2.m');figure(1);cho_display(net);

uses mumpsx.net

a=pi/4

anchor p1 [A]   [l=10u w=10u oz=pi]
beam3d p1 [A B] [l=200u w=2u h=2u oz=0]

beamjoint12 p1 [B C] [l=200u w=2u h=2u oz=a oy=0.011]

beam3d p1 [C D] [l=200u w=2u h=2u oz=0]
anchor p1 [D]   [l=10u w=10u oz=pi]

f3d * [B][F=10e-6 oy=pi/2] 

