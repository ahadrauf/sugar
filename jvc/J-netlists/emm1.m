% net=cho_load('emm1.m');    figure(1);   cho_display(net); shownodes(net);
% q(lookup_coord(net,'a2','rz'))
%net=cho_load('emm1.m',param);   [q]=cho_dc(net);    figure(1);   cho_display(net,q*0);

uses soi.net
param H = 2u
param L = 100u
param W = 2u
param L_anchor = 100u

%Middle anchor and beams
anchor p2 [a1]    [l=10u w=10u h=H oz=pi/4]
beam3d p1 [a1 a2] [l=L   w=W   h=H oz=pi Youngsmodulus=100e9]
beam3d p1 [a1 b1] [l=L_anchor  w=W   h=H oz=-pi/2]
beam3d p1 [b1 b2] [l=L   w=W   h=H oz=pi]
anchor p1 [b1]    [l=10u w=10u h=H oz=0]

f3d * [a2][F=-600u oz=pi/2] 
