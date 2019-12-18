% net=cho_load('can0.m');q=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','y'))

uses mumps2.net
uses process_cmos1.m

%SG=4000
anchor  p1 [A]   [l=20u  w=20u      oz=pi]
beam3d_a1  p1 [A B] [l=100u w=2u  h=2u         straingradient=0]
%beam3d  p1 [B C] [l=100u w=2u  h=2u oz=pi/2 straingradient=SG]
%beam3d  p1 [C D] [l=100u w=2u  h=2u oz=pi   straingradient=SG]
%anchor  p1 [C]     [l=10u w=10u]
f3d * [B] [F=4u oz=pi/2] 

%resistor * [a b] [l=100u]

