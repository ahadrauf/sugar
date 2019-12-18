%beam and slider 
% net=cho_load('_circ3.m');[q]=cho_dc(net);figure(1);cho_display(net,q); 
% net=cho_load('_circ3.m');figure(1);cho_display(net);  
 
uses mumps.net 
 
pi=3.1415926535 
%anchor   p1 [a]  [l=1u w=1u h=1u oz=-90] 

a=pi/2
a1=pi/3
w=10u
h=10u
R1=20u
Rc=100u
R2=Rc*2
wc=Rc*1.5
Rcb=wc/2
z1=0
aL=50u
aW=50u
aH=h

%1st weave
z=z1
circbeam p1 [a02 a03][w=w h=h radius=R1 alpha=a oz=z ox=180] %small half
circbeam p1 [a03 a04][w=w h=h radius=R2+R1 alpha=a/2 oz=z-90 ox=180] %half curve
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
anchor   p1 [a17]    [l=aL w=aW h=aH oz=z] 

%2nd weave
z=z1+120
circbeam p1 [b02 b03][w=w h=h radius=R1 alpha=a oz=z ox=180] %small half
circbeam p1 [b03 b04][w=w h=h radius=R2+R1 alpha=a/2 oz=z-90 ox=180] %half curve
circbeam p1 [b04 b05][w=w h=h radius=R1 alpha=a oz=z-90-45 ox=0] %small loop up
circbeam p1 [b05 b06][w=w h=h radius=R1 alpha=a oz=z-45 ox=0] 
circbeam p1 [b06 b07][w=w h=h radius=R2+2*R1 alpha=a oz=z-45+90 ox=0] %curve A
circbeam p1 [b07 b08][w=w h=h radius=R1 alpha=a oz=z+180-45 ox=180] %small loop dn
circbeam p1 [b08 b09][w=w h=h radius=R1 alpha=a oz=z+90-45 ox=180] 
circbeam p1 [b09 b10][w=w h=h radius=R2+4*R1 alpha=a oz=z-45 ox=180] %curve B
circbeam p1 [b10 b11][w=w h=h radius=R1 alpha=a oz=z-90-45 ox=0] %small loop up
circbeam p1 [b11 b12][w=w h=h radius=R1 alpha=a oz=z-45 ox=0] 
circbeam p1 [b12 b13][w=w h=h radius=R2+6*R1 alpha=a oz=z-45+90 ox=0] %curve C
circbeam p1 [b13 b14][w=w h=h radius=R1 alpha=a oz=z+180-45 ox=180] %small loop dn
circbeam p1 [b14 b15][w=w h=h radius=R1 alpha=a oz=z+90-45 ox=180] 
circbeam p1 [b15 b16][w=w h=h radius=R2+8*R1 alpha=a/2 oz=z-45 ox=180] %curve D
circbeam p1 [b16 b17][w=w h=h radius=R1 alpha=a oz=z-90 ox=0] %small half
anchor   p1 [b17]    [l=aL w=aW h=aH oz=z] 

%3rd weave
z=z1+240
circbeam p1 [c02 c03][w=w h=h radius=R1 alpha=a oz=z ox=180] %small half
circbeam p1 [c03 c04][w=w h=h radius=R2+R1 alpha=a/2 oz=z-90 ox=180] %half curve
circbeam p1 [c04 c05][w=w h=h radius=R1 alpha=a oz=z-90-45 ox=0] %small loop up
circbeam p1 [c05 c06][w=w h=h radius=R1 alpha=a oz=z-45 ox=0] 
circbeam p1 [c06 c07][w=w h=h radius=R2+2*R1 alpha=a oz=z-45+90 ox=0] %curve A
circbeam p1 [c07 c08][w=w h=h radius=R1 alpha=a oz=z+180-45 ox=180] %small loop dn
circbeam p1 [c08 c09][w=w h=h radius=R1 alpha=a oz=z+90-45 ox=180] 
circbeam p1 [c09 c10][w=w h=h radius=R2+4*R1 alpha=a oz=z-45 ox=180] %curve B
circbeam p1 [c10 c11][w=w h=h radius=R1 alpha=a oz=z-90-45 ox=0] %small loop up
circbeam p1 [c11 c12][w=w h=h radius=R1 alpha=a oz=z-45 ox=0] 
circbeam p1 [c12 c13][w=w h=h radius=R2+6*R1 alpha=a oz=z-45+90 ox=0] %curve C
circbeam p1 [c13 c14][w=w h=h radius=R1 alpha=a oz=z+180-45 ox=180] %small loop dn
circbeam p1 [c14 c15][w=w h=h radius=R1 alpha=a oz=z+90-45 ox=180] 
circbeam p1 [c15 c16][w=w h=h radius=R2+8*R1 alpha=a/2 oz=z-45 ox=180] %curve D
circbeam p1 [c16 c17][w=w h=h radius=R1 alpha=a oz=z-90 ox=0] %small half
anchor   p1 [c17]    [l=aL w=aW h=aH oz=z] 

%center
z=z1
beam3d   p1 [a01 a02][l=Rcb w=w h=h oz=z ox=0] 
beam3d   p1 [b01 b02][l=Rcb w=w h=h oz=z+120 ox=0] 
beam3d   p1 [c01 c02][l=Rcb w=w h=h oz=z+240 ox=0] 
circbeam p1 [a01 b01m][w=wc h=h radius=Rc alpha=a1 oz=z+90 ox=0] 
circbeam p1 [b01m b01][w=wc h=h radius=Rc alpha=a1 oz=z+90+60 ox=0] 
circbeam p1 [b01 c01m][w=wc h=h radius=Rc alpha=a1 oz=z+90+60*2 ox=0] 
circbeam p1 [c01m c01][w=wc h=h radius=Rc alpha=a1 oz=z+90+60*3 ox=0] 
circbeam p1 [c01 a01m][w=wc h=h radius=Rc alpha=a1 oz=z+90+60*4 ox=0] 
circbeam p1 [a01m a01][w=wc h=h radius=Rc alpha=a1 oz=z+90+60*5 ox=0] 

f=200u
f3d * [a01][F=f oy=90] 
f3d * [b01][F=f oy=90] 
f3d * [c01][F=f oy=90] 
  
