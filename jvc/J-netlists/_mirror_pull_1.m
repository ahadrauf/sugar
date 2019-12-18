%beam and slider
% net=cho_load('_mirror_pull_1.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('_mirror_pull_1.m');figure(1);cho_display(net);
% net=cho_load('_mirror_pull_1.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);
% net=cho_load('_mirror_pull_1.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);q(lookup_coord(net,'g','ry'))
% q=zeros(12,1);q(4)=pi/4;figure(1);cho_display(net,q)

uses mumps.net
W=40u
H=10u
anchor p1 [a][l=20u w=20u h=20u oz=-pi/2]
%beam3d p1 [c1 a][l=600u w=20u h=2u oz=0]
%beam3d p1 [c1 b][l=600u w=20u h=2u oz=pi]

Wcirc=3u
Hcirc=2u
OZcirc=pi/4
Rcirc=200u
circbeam3 p1 [a c2][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=0]
%circbeam3 p1 [a c][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=0]
%beam3d p1 [a b][l=Rcirc*OZcirc w=Wcirc h=Hcirc oz=0]
%circbeam3 p1 [c2 c3][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=pi/2+OZcirc*1]
%circbeam3 p1 [c3 c4][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=pi/2+OZcirc*2]
%circbeam3 p1 [c4 c5][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=pi/2+OZcirc*3]
%circbeam3 p1 [c5 c6][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=pi/2+OZcirc*4]
%circbeam3 p1 [c6 c7][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=pi/2+OZcirc*5]
%circbeam3 p1 [c7 c8][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=pi/2+OZcirc*6]
%circbeam3 p1 [c8 c9][w=Wcirc h=Hcirc radius=Rcirc alpha=OZcirc oz=pi/2+OZcirc*7]

%circbeam3 p1 [d dd][w=Wcirc h=Hcirc radius=100u alpha=pi/4 oz=pi/2]

f3d * [c2][M=1e-7 oy=pi/2]
%f3d * [b][F=-40u]
%rigidlinkbeamcorner p1 [B C] [l=100u w=10u L1=10u oz1=0 L2=10u oz2=0 oz=pi/4*0]
