% net=cho_load('nox1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'D','y'))
% net=cho_load('nox1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('nox1.m');figure(1);cho_display(net);

uses mumpsx.net
anchor p1 [A]   [l=10u w=10u oz=180]
beam3d p1 [A B] [l=100u w=2u h=2u]
beam3d p1 [B C] [l=100u w=2u h=2u oz=90]
beam3d p1 [C D] [l=100u w=2u h=2u]
f3d * [D][F=1e-6 oz=90] 

%beam3d p1 [A E] [l=100u w=2u h=2u]
%beam3d p1 [E DD] [l=100u w=2u h=2u]
%f3d * [E][M=(1e-6 * 100e-6) oy=90] 
%f3d * [DD][F=1u oz=90] 

