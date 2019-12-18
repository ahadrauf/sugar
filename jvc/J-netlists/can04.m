%beam and slider
% net=cho_load('can04.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('can04.m');figure(1);cho_display(net);

uses mumps.net

a=20
b=90
d=50u
anchor p1 [e]   [l=1u w=1u h=1u oz=-90]
beam3d p1 [a b][l=180u w=d*2 h=10u oy=-a]
slider2 p1 [b c][l=120u w=d*2 h=10u oy=-a slideroz=0]

beam3d p1 [a al][l=d w=1u h=1u oz=-90]
beam3d p1 [a ar][l=d w=1u h=1u oz=90]
beam3d p1 [ar aar][l=100u w=2u h=2u oy=-90]
beam3d p1 [al aal][l=100u w=2u h=2u oy=-90]
beam3d p1 [e aal][l=d w=1u h=1u oz=-90]
beam3d p1 [e aar][l=d w=1u h=1u oz=90]
beam3d p1 [e1 aal][l=d w=10u h=3u oz=-90 ox=90]
beam3d p1 [e1 aar][l=d w=10u h=3u oz=90 ox=90]

m=1e-9
n=-5u
%n=-0.001u
f3d * [a][F=n oz=0] 
f3d * [c][F=n oz=0] 
%f3d * [ar][F=n oz=0] 
%f3d * [a][F=-300u oy=90] 
%f3d * [al][M=-m oz=90] 
%f3d * [ar][M=-m oz=90] 
%f3d * [aal][M=-m oz=90] 
%f3d * [aar][M=m oz=-90] 

c=160
beam3d p1 [a cal][l=d w=1u h=1u oz=-90]
beam3d p1 [a car][l=d w=1u h=1u oz=90]
beam3d p1 [car caar][l=150u w=2u h=2u oy=-c]
beam3d p1 [cal caal][l=150u w=2u h=2u oy=-c]
beam3d p1 [ce caal][l=d w=1u h=1u oz=-90]
beam3d p1 [ce caar][l=d w=1u h=1u oz=90]
beam3d p1 [ce1 caal][l=d w=10u h=3u oz=-90 ox=-c]
beam3d p1 [ce1 caar][l=d w=10u h=3u oz=90 ox=c]

%a=60
%b=a*3.1415/180
%anchor p1 [C]   [l=10u w=10u oz=-90]
%beam3d p1 [A C] [l=200u w=2u h=2u oz=-a]
%slider2 p1 [A B] [l=200u w=2u h=2u oz=a slideroz=b]
%beam3d p1 [B BB] [l=10u w=10u h=10u oz=a]

%f3d * [A][F=-30u] 

%beam3d p1 [A1 C]  [l=200u w=1e-22 h=1u oz=-a]
%beam3d p1 [A1 B1] [l=50u w=1e-22 h=1u oz=a]
%beam3d p1 [B1 B2] [l=200u w=4u h=1u oz=a]
