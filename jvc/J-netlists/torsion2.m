% net=cho_load('torsion2.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('torsion2.m');[q]=cho_dc(net);figure(1);cho_display(net,q);qy=q(lookup_coord(net,'B','y')),qz=q(lookup_coord(net,'B','z'))
% net=cho_load('torsion2.m');figure(1);cho_display(net);

uses mumpsx.net
anchor p1 [A]   [l=10u w=10u oz=180]
beam3d p1 [A Bl] [l=50u w=19u h=2u]
beam3d p1 [Bl B] [l=50u w=19u h=2u ox=90]
beam3d p1 [B Br] [l=50u w=19u h=2u ox=90]
beam3d p1 [Br C] [l=50u w=19u h=2u]
anchor p1 [C]   [l=10u w=10u oz=0]

%f3d * [B][F=50000u oy=90] 
f3d * [B][F=50000u oz=90] 
%f3d * [B][M=3 ] 

