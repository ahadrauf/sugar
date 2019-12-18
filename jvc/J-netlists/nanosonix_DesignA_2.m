% net=cho_load('nanosonix_DesignA_2.m',p); q=cho_dc(net); figure(1); clf; cho_display(net,q); 
%one = q(lookup_coord(net, 'c', 'y'))

uses mumps.net
L = 100u
w = 10u
L1 = 4u
hp1 = 2u
ho2 = 2u
hp0 = 0.5u
Lanchor = 15u
wanchor = Lanchor

%anchor / oxide
anchor p1 [a] [l=Lanchor w=wanchor oz=-pi]
beam3dlink p1 [a c] [L1=Lanchor/2 l=(hp1/2 + ho2 + hp0) w=wanchor h=wanchor oy1=pi/2 oy=pi/2]

%Cantilever
beam3d p1 [a b] [l=L w=w ]
f3d * [b] [ M=1p oz=pi/2 ]

%Electrodes
beam3dlink p0 [c d] [L1=L1+Lanchor/2 l=L/4-L1 w=w*3 ]
beam3dlink p0 [d e] [L1=L1 l=2*L/4-L1 w=w*3 ]
beam3dlink p0 [e f] [L1=L1 l=L/4-L1 w=w*3 L2=L1]

%Tip electrode
beam3dlink p1 [f g] [L1=(hp1/2 + ho2 + hp0/2) oy1=-pi/2 l=L/10 w=w ]
anchor p1 [g] [l=15u w=15u ]
beam3dlink p1 [g h] [L1=Lanchor/2 l=(hp1/2 + ho2 + hp0) w=wanchor h=wanchor oy1=-pi/2 oy=pi/2]


