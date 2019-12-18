%net=cho_load('perforated3.m');figure(1);cho_display(net);
%net=cho_load('perforated3.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net
uses perforatedcomb4.m
%anchor p1 [left][l=10u w=10u h=10u oz=-pi]
%perforatedcomb p1 [left right][numberfingers=4 w=20u h=10u l=50u whorizontal=2u wfinger=2u lfinger=40u gap=2u] 

%f3d * [right][M=0.05u]
%f3d * [right][F=-1000u oy=pi/2]

perforatedcomb4 p1 [n1 n2][numberfingers=4 w=20u h=10u whorizontal=6u wfinger=3u lfinger=30u gap=3u wvertical=6u L1=20u L2=20u] 
beam3d p1 [n1 nn][l=50u w=10u oz=pi/4]
beam3d p1 [n2 nnn][l=50u w=10u oz=pi/4]
