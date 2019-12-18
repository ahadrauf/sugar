%beam and slider
% net=cho_load('_circ10.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('_circ10.m');figure(1);cho_display(net);
% net=cho_load('_circ10.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('_circ10.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);
% net=cho_load('_circ10.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);q(lookup_coord(net,'g','ry'))
% q=zeros(12,1);q(4)=pi/4;figure(1);cho_display(net,q)

uses mumps.net
W=10u
H=2u
anchor   p1 [d]  [l=10u w=2u h=2u oz=-pi/2]
%anchor   p1 [M]  [l=10u w=2u h=2u oz=-pi/2]
%beam3d   p1 [d M] [l=20u w=2u h=2u oy=-pi/2]

circbeam10 p1 [b d][w=W h=H radius=100u alpha=pi/3 oz=-pi/3]
circbeam10 p1 [d a][w=W h=H radius=100u alpha=pi/3 oz=0]

theta=(pi/4)/3
R=100u 
beam3d p1 [d  a1] [l=R*theta w=W h=H oz=(1-cos(theta))/theta]
beam3d p1 [a1 a2] [l=R*theta w=W h=H oz=(1-cos(theta*3))/theta/3*1.05]
beam3d p1 [a2  c] [l=R*theta w=W h=H oz=(1-cos(theta*9))/theta/9*0.9]

oz=pi
beam3d p1 [d  b1] [l=R*theta w=W h=H oz=oz-(1-cos(theta))/theta]
beam3d p1 [b1 b2] [l=R*theta w=W h=H oz=oz-(1-cos(theta*3))/theta/3*1.05]
beam3d p1 [b2 e] [l=R*theta w=W h=H oz=oz-(1-cos(theta*9))/theta/9*0.9]

oz=pi/2
f3d * [a] [F=10000u oz=oz]
f3d * [b] [F=10000u oz=oz]
f3d * [c] [F=10000u oz=oz]
f3d * [e] [F=10000u oz=oz]

%f3d * [d] [F=500u oy=pi/2]
%f3d * [M] [F=500u oy=pi/2]

%f3d * [d] [M=0.010000u oz=0*pi/2]
%f3d * [M] [M=0.010000u oz=0*pi/2]

