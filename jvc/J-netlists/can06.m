% p.fx=fx;p.fy=fy;p.m=m; net=cho_load('can06.m',p); q=cho_dc(net); figure(1); cho_display(net,q); qx=q(lookup_coord(net,'B','x')),qy=q(lookup_coord(net,'B','y')),qrz=q(lookup_coord(net,'B','rz'))

uses mumps2.net

param fx
param fy
param m

anchor  p1 [A]   [l=20u  w=20u oz=pi]
beam3d  p1 [A B] [l=50u w=2u  h=2u ]

f3d * [B] [F=fx oz=0] 
f3d * [B] [F=fy oz=pi/2] 
f3d * [B] [M=m oy=-pi/2] 

