%net=cho_load('comb1.m');figure(1);cho_display(net);
%net=cho_load('comb1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net
uses comb.m
anchor p1 [middle][l=10u w=10u h=10u oz=pi/2]
comb p1 [left middle right][V=0 numberfingers=11 gap=9u fingerwidth=3u fingerheight=50u fingerlength=30u supportwidth=3u supportheight=50u oz=pi/2*0] 


