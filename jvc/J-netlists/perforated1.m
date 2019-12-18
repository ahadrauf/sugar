%net=cho_load('perforated1.m');figure(1);cho_display(net);
%net=cho_load('perforated1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net
uses perforated.m
anchor p1 [left][l=10u w=10u h=10u oz=pi/2]
perforated p1 [left right][numberhole=5 w=20u h=20u l=50u whorizontal=2u wvertical=2u] 


