% net=cho_load('predisp1A.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% p.qx=q(1);p.qy=q(2);p.qz=q(3);p.qox=q(4);p.qoy=q(5);p.qoz=q(6);net=cho_load('predisp1A.m',p);[q]=cho_dc(net);figure(1);cho_display(net,q);
% p.qx=q(1);p.qy=q(2);p.qz=q(3);p.qox=q(4);p.qoy=q(5);p.qoz=q(6);net=cho_load('predisp1A.m',p);[qq]=cho_dc(net);figure(1);cho_display(net,qq);

uses mumpsx.net
uses predisplacedbeam3a.net
param qx=0 
param qy=110u 
param qz=0 
param qox=0 
param qoy=0 
param qoz=0

anchor  p1 [B]   [l=10u w=10u h=10u oz=-pi/2]

%beam3d  p1 [B BB] [l=200u w=5u h=2u oz=0]
%f3d     *  [BB]   [F=60u oz=90]

A=1
B=1
C=1

if A==B
    [
        predisplacedbeam3a p1 [B C][l=200u w=5u oz=0 qx=qx qy=qy qz=qz qox=qox qoy=qoy qoz=qoz]
    ]
if A==C
    [
        predisplacedbeam3a p1 [B C][l=200u w=50u oz=0 qx=qx qy=qy qz=qz qox=qox qoy=qoy qoz=qoz]
    ]
%prestressedbeam   p1 [B C][l=200u w=5u oz=0 qx=qx qy=qy qz=qz qox=qox qoy=qoy qoz=qoz]

beam3d  p1 [C d] [l=100u w=2u h=2u oz=0]
beam3d  p1 [d e] [l=50u w=2u h=2u oz=-pi/2]
anchor  p1 [e]   [l=10u w=10u h=10u oz=-pi/2]



%beam3d  p1 [B CC] [l=200u w=5u h=2u oz=0]
%prestressedbeam  p1 [B CC][l=200u w=5u oz=0 qx=qx qy=qy qz=qz qox=qox*180/pi qoy=qoy*180/pi qoz=qoz*180/pi]


%predisplacedbeam3 p1 [B CCC][l=200u w=5u oz=0 qx=0*qx qy=0*qy qz=0*qz qox=0*qox*180/pi qoy=0*qoy*180/pi qoz=0*qoz*180/pi]
%prestressedbeam   p1 [B CCC][l=200u w=5u oz=0 qx=-qx qy=-qy qz=-qz qox=-qox*180/pi qoy=-qoy*180/pi qoz=-qoz*180/pi]


