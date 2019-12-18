%net=cho_load('test2.m');figure(1);cho_display(net);
%net=cho_load('test2.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net
uses subnet_comb.m
anchor p1 [middle][l=10u w=10u h=10u oz=pi/2]
subnet_comb p1 [left middle right][V=1000 numberfingers=11 gap=6u fingerwidth=2u fingerheight=2u fingerlength=20u supportwidth=2u supportheight=2u oz=pi/2*0] 


