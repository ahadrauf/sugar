% net=cho_load('_circ12.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

uses mumps.net
anchor   p1 [a]  [l=0.001u w=0.001u h=0.001u oz=-pi/2]

W = 5u
H = 2u
R = 100u
circbeam12 p1 [a b][w=W h=H radius=R alpha=-pi/2 oz=-pi/2 Youngsmodulus=170e9] 
circbeam12 p1 [a bx][w=W h=H radius=R alpha=-pi/2 oz=-pi/2] 

oz=0
ox=0
oy=pi/2
f=7u

f3d * [b]  [F=f oz=oz ox=ox oy=oy] 

