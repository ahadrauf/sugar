% p.theta=10; net=cho_load('ETA1.m',p); q=cho_dc(net); figure(1); clf; cho_display(net,q); one = q(lookup_coord(net, 'c', 'y'))


uses mumps.net
param L = 200u
w = 2u
param theta = 10
theta = theta * pi/180
T = 310

anchor p1 [a] [l=10u w=10u oz=-pi/2]
beam3d p1 [a c] [l=L w=w oz=theta T=T]
beam3d p1 [c b] [l=L w=w oz=-theta T=T]
anchor p1 [b] [l=10u w=10u oz=-pi/2]


%beam3d p1 [c d] [l=L w=w oz=pi/2-theta ]
%beam3d p1 [d e] [l=L w=w oz=pi/2+theta ]
%anchor p1 [e] [l=10u w=10u oz=pi/2]

