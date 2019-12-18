% net=cho_load('mass1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('mass1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('mass1.m');figure(1);cho_display(net);

uses mumps.net
mL = 500u
H = 20u
W = 2u
W2 = 20u
fL = 200u
aL = 50u 

beam2d p1 [c d][l=mL w=mL h=H oz=0] 

anchor p1 [a][l=10u w=10u oz=pi/2]
beam2d p1 [a c][l=fL w=W h=H oz=-pi/2] 

beam2d p1 [c e][l=fL w=W h=H oz=-pi/2] 
anchor p1 [e][l=10u w=10u oz=-pi/2]

f3d * [c][F=100e-6 oz=0] 
