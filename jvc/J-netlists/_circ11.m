%beam and slider
% net=cho_load('_circ11.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('_circ11.m');figure(1);cho_display(net);
% net=cho_load('_circ11.m');[q,k]=cho_dc3(net);figure(1);cho_display(net,q);q(lookup_coord(net,'g','ry'))
% q=zeros(12,1);q(4)=pi/4;figure(1);cho_display(net,q)
% net=cho_load('_circ11.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net
W=10u
H=2u
anchor   p1 [d]  [l=10u w=2u h=2u oz=-pi/2]
%anchor   p1 [M]  [l=10u w=2u h=2u oz=-pi/2]
%beam3d   p1 [d M] [l=20u w=2u h=2u oy=-pi/2]

%circbeam11 p1 [b d][w=W h=H radius=100u alpha=pi/4 oz=-pi/4] 
%circbeam11 p1 [d a][w=W h=H radius=100u alpha=pi/4 oz=0] 

circbeam12 p1 [d b2][w=W h=H radius=100u alpha=pi/4 ] 
%circbeam12 p1 [b2 d] [w=W h=H radius=100u alpha=pi/4/2 oz=-pi/8] 

%circbeam12 p1 [d a1] [w=W h=H radius=100u alpha=pi/4/2 oz=0] 
%circbeam12 p1 [a1 a2][w=W h=H radius=100u alpha=pi/4/2 oz=pi/8] 

theta=(pi/4)/3
R=100u 
%beam3d p1 [d  a1] [l=R*theta w=W h=H oz=(1-cos(theta))/theta] beam3d p1 [a1 a2] [l=R*theta w=W h=H oz=(1-cos(theta*3))/theta/3*1.05] beam3d p1 [a2  c] [l=R*theta w=W h=H oz=(1-cos(theta*9))/theta/9*0.9]
oz=pi
%beam3d p1 [d  b1] [l=R*theta w=W h=H oz=oz-(1-cos(theta))/theta] beam3d p1 [b1 b2] [l=R*theta w=W h=H oz=oz-(1-cos(theta*3))/theta/3*1.05] beam3d p1 [b2 e] [l=R*theta w=W h=H oz=oz-(1-cos(theta*9))/theta/9*0.9]

oz=0*pi/2
ox=0
oy=-pi/2
f=1000u
%f3d * [b]  [F=f oz=oz ox=ox oy=oy] 
%f3d * [a]  [F=f oz=oz ox=ox oy=oy]
%f3d * [a2] [F=f oz=oz ox=ox oy=oy] 
%f3d * [b1] [F=f oz=oz ox=ox oy=oy] 

%f3d * [e] [F=f oz=oz ox=ox oy=oy] f3d * [c] [F=f oz=oz ox=ox oy=oy]
