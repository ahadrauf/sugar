%net=cho_load('perforated2.m');figure(1);cho_display(net);
%net=cho_load('perforated2.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net
uses perforated.m
%anchor p1 [left][l=10u w=10u h=10u oz=-pi]
perforated p1 [left right][numberholes=4 w=20u h=10u l=50u whorizontal=2u wvertical=2u] 

%f3d * [right][M=0.05u]
%f3d * [right][F=-1000u oy=pi/2]
