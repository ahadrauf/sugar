% net=cho_load('distload1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('distload1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('distload1.m');figure(1);cho_display(net);

uses mumps.net
anchor p1 [A]   [l=10u w=10u oz=180]
beam3d p1 [A B] [l=200u w=5u h=20u]
beam3d p1 [B C] [l=200u w=5u h=20u]
p3d * [A B][l=100u P="70"] 
p3d * [B C][l=100u P="-250000*x"] 
