%clear all; i=0; p.T=310; p.L=200e-6; for theta = 0:1/20:10, i=i+1; p.theta=theta; net=cho_load('ETA2.m',p);  q=cho_dc(net); Y(i) = q(lookup_coord(net, 'd','y')); X(i) = theta; end, figure(1); hold on; grid on; plot(X,Y,'b');

uses mumps.net
param L = 200u
param w = 2u
param T = 310
param a = 0.5 %Ratio (L*a)/(L*(1-a))
param theta = 10
theta = theta * pi/180

anchor p1 [a] [l=10u w=10u oz=pi]
beam3d p1 [a b] [l=L*a w=w T=T]
beam3d p1 [b c] [l=L*atan(theta) w=w T=T oz=pi/2]
beam3d p1 [c d] [l=L*(1-a) w=w T=T]

beam3d p1 [d u] [l=10u w=w*4 oz=pi/2]
beam3d p1 [d n] [l=10u w=w*4 oz=-pi/2]

beam3d p1 [d e] [l=L*(1-a) w=w T=T]
beam3d p1 [e f] [l=L*atan(-theta) w=w T=T oz=pi/2]
beam3d p1 [f g] [l=L*a w=w T=T]
anchor p1 [g] [l=10u w=10u ]


