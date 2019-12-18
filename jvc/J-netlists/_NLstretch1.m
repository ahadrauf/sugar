%beam and slider
% p.f=600e-6;net=cho_load('_NLstretch1.m',p);q=cho_dc_secant(net);figure(1);cho_display(net,q);q(2),q(2+6)

uses mumps.net

param f=100u

anchor   p1 [a]  [l=10u w=10u h=10u oz=pi]

%nls_beam p1 [a b][w=2u h=2u l=50u]
%nls_beam p1 [b c][w=2u h=2u l=50u]

beam   p1 [a B][w=2u h=2u l=50u]
beam   p1 [B c][w=2u h=2u l=50u]

anchor   p1 [c]  [l=10u w=10u h=10u ]

%f3d * [b][F=f oz=-pi/2]
f3d * [B][F=f oz=-pi/2]

