% net=cho_load('predispP1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% p.qx=q(1);p.qy=q(2);p.qz=q(3);p.qox=q(4);p.qoy=q(5);p.qoz=q(6);net=cho_load('predisp1.net',p);[q]=cho_dc(net);figure(1);cho_display(net,q);
%p.qx=q(1);p.qy=q(2);p.qz=q(3);p.qox=q(4);p.qoy=q(5);p.qoz=q(6);net=cho_load('predisp1.net',p);[qq]=cho_dc(net);figure(1);cho_display(net,qq);

%predisp_unimorph1
%net=cho_load('predisp_unimorph1.m');[qq]=cho_dc(net);figure(1);cho_display(net,qq);dX=qq(lookup_coord(net,'m1','x')), dY=qq(lookup_coord(net,'m1','y')), dZ=qq(lookup_coord(net,'m1','z')), RX=qq(lookup_coord(net,'m1','rx')), RY=qq(lookup_coord(net,'m1','ry')), RZ=qq(lookup_coord(net,'m1','rz'))

uses mumpsx.net
uses predisplacedbeam3a.net
pi=3.141592
param qx=0 
param qy=110u 
param qz=0 
param qox=0 
param qoy=0 
param qoz=0
W = 200u

anchor  p1 [A]    [l=10u w=10u h=10u oz=pi]
%right side
beam3d  p1 [A AR] [l=200u w=W h=2u oz=0]
predisplacedbeam3a p1 [AR m1][l=100u w=W oz=0   qx=0 qy=0 qz=2u   qox=0 qoy=-pi/2*1.75 qoz=0]
%left side
beam3d  p1 [A AL] [l=200u w=W h=2u oz=pi]
predisplacedbeam3a p1 [AL m2][l=100u w=W oz=pi   qx=0 qy=0 qz=2u   qox=0 qoy=-pi/2*1.75 qoz=0]

f3d     *  [m1]   [F=600u oy=-pi/2]
%f3d     *  [m1]   [M=0.001u oz=0]


%prestressedbeam   p1 [B C][l=200u w=5u  oz=0 qx=0 qy=1u qz=0 qox=0 qoy=0 qoz=45*pi/180]

%predisplacedbeam3a p1 [C D][l=200u w=5u oz=45 qx=0 qy=1u qz=0 qox=0 qoy=0 qoz=45]
%prestressedbeam    p1 [C D][l=200u w=5u oz=45 qx=0 qy=1u qz=0 qox=0 qoy=0 qoz=45*pi/180]
%f3d     *  [C]   [F=-20000u oz=0]
%prestressedbeam   p1 [B C][l=200u w=5u oz=0 qx=qx qy=qy qz=qz qox=qox*180/pi qoy=qoy*180/pi qoz=qoz*180/pi]
%predisplacedbeam3a p1 [B C][l=200u w=5u oz=0 qx=qx qy=qy qz=qz qox=qox*180/pi qoy=qoy*180/pi qoz=qoz*180/pi]
%beam3d  p1 [C d] [l=100u w=2u h=2u oz=0]
%beam3d  p1 [d e] [l=50u w=2u h=2u oz=-90]
%anchor  p1 [e]   [l=10u w=10u h=10u oz=-90]
%beam3d  p1 [C d] [l=47u w=80u h=2u oz=-90 ox=90]
%beam3d  p1 [c C] [l=47u w=80u h=2u oz=-90 ox=90]
%f3d     *  [c]   [F=-2000u oz=0]
%anchor  p1 [d]   [l=10u w=10u h=10u oz=-90]
%beam3d  p1 [B CC] [l=200u w=5u h=2u oz=0]
%prestressedbeam  p1 [B CC][l=200u w=5u oz=0 qx=qx qy=qy qz=qz qox=qox*180/pi qoy=qoy*180/pi qoz=qoz*180/pi]
%predisplacedbeam3 p1 [B CCC][l=200u w=5u oz=0 qx=0*qx qy=0*qy qz=0*qz qox=0*qox*180/pi qoy=0*qoy*180/pi qoz=0*qoz*180/pi]
%prestressedbeam   p1 [B CCC][l=200u w=5u oz=0 qx=-qx qy=-qy qz=-qz qox=-qox*180/pi qoy=-qoy*180/pi qoz=-qoz*180/pi]


