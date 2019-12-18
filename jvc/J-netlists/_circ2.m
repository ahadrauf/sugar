%beam and slider 
% net=cho_load('_circ2.m');[q]=cho_dc(net);figure(1);cho_display(net,q); 
% net=cho_load('_circ2.m');figure(1);cho_display(net);  
 
uses mumps.net 
 
pi=3.1415926535 
%anchor   p1 [a]  [l=1u w=1u h=1u oz=-90] 

a=pi/2
a1=pi/3
w=10u
h=10u
R2=100u
R1=20u
Rc=100u
wc=Rc
%first weave
z=45
circbeam p1 [a02 a03][w=w h=h radius=R1 alpha=a oz=z ox=180] %small half
circbeam p1 [a03 a04][w=w h=h radius=R2 alpha=a/2 oz=z-90 ox=180] %half curve
circbeam p1 [a04 a05][w=w h=h radius=R1 alpha=a oz=z-90-45 ox=0] %small loop up
circbeam p1 [a05 a06][w=w h=h radius=R1 alpha=a oz=z-45 ox=0] 
circbeam p1 [a06 a07][w=w h=h radius=R2+2*R1 alpha=a oz=z-45+90 ox=0] %curve A
circbeam p1 [a07 a08][w=w h=h radius=R1 alpha=a oz=z+180-45 ox=180] %small loop dn
circbeam p1 [a08 a09][w=w h=h radius=R1 alpha=a oz=z+90-45 ox=180] 
circbeam p1 [a09 a10][w=w h=h radius=R2+4*R1 alpha=a oz=z-45 ox=180] %curve B
circbeam p1 [a10 a11][w=w h=h radius=R1 alpha=a oz=z-90-45 ox=0] %small loop up
circbeam p1 [a11 a12][w=w h=h radius=R1 alpha=a oz=z-45 ox=0] 
circbeam p1 [a12 a13][w=w h=h radius=R2+6*R1 alpha=a oz=z-45+90 ox=0] %curve C
circbeam p1 [a13 a14][w=w h=h radius=R1 alpha=a oz=z+180-45 ox=180] %small loop dn
circbeam p1 [a14 a15][w=w h=h radius=R1 alpha=a oz=z+90-45 ox=180] 
circbeam p1 [a15 a16][w=w h=h radius=R2+8*R1 alpha=a/2 oz=z-45 ox=180] %curve D
circbeam p1 [a16 a17][w=w h=h radius=R1 alpha=a oz=z-90 ox=0] %small half
anchor   p1 [a17]    [l=20u w=20u h=20u oz=z] 

beam3d   p1 [a01 a02][l=Rc/2 w=w h=h oz=z ox=0] 
circbeam p1 [a01 b01m][w=wc h=h radius=Rc alpha=a1 oz=z+90 ox=0] 
circbeam p1 [b01m b01][w=wc h=h radius=Rc alpha=a1 oz=z+90+60 ox=0] 
circbeam p1 [b01 c01m][w=wc h=h radius=Rc alpha=a1 oz=z+90+60*2 ox=0] 
circbeam p1 [c01m c01][w=wc h=h radius=Rc alpha=a1 oz=z+90+60*3 ox=0] 
circbeam p1 [c01 a01m][w=wc h=h radius=Rc alpha=a1 oz=z+90+60*4 ox=0] 
circbeam p1 [a01m a01][w=wc h=h radius=Rc alpha=a1 oz=z+90+60*5 ox=0] 

f3d * [a02][F=1 ] 
  
