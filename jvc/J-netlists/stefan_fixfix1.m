%  net = cho_load('stefan_fixfix1.m');dq = cho_dc(net); defl = dqval(net,dq,'b','y')

uses mumps.net
uses stdlib.net

anchor  p1 [a]   [l=5u w=4u oz=deg(-90)]
beam2d  p1 [a b] [l=50u w=1u]
beam2d  p1 [b c] [l=50u w=1u]
anchor  p1 [c]   [l=5u w=4u oz=deg(-90)]
f3d     *  [b]   [F=5u oz=deg(-90)]

