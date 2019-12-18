% net=cho_load('beamplate.m');q=cho_dc(net);figure(1);cho_display(net,q); 

uses mumps2.net

anchor p1 [UL][l=10u w=10u h=10u oz=180] 

beam3da p1 [UL B][l=200u w=2u h=20u oz=0]

f3d * [B][ F=10u oz=90] 

