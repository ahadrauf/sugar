%net=cho_load('jdrive_mg_1.m');q=cho_dc(net);figure(1);cho_display(net,q);y=q(18-6+2)
%uses memgen.m

uses mumps2.net

H = 24u
space = 27u

anchor p1   [a]     [l=100u      w=100u h=240u oz=-pi] 

beam3d p1  [a b]   [l=1200u    w=24u h=H oz=pi/2]

beam3d p1  [b c]   [l=space      w=24u h=H oz=0]

beam3d p1  [c d]   [l=1200u    w=24u h=H oz=-pi/2]

beam3d p1  [d e]   [l=1200u      w=24u h=240u oz=0] %actuator

anchor p1   [e]     [l=100u      w=100u h=240u oz=0] 

V=7000
eps=8.854e-12
L = 1200u
W = 200u
area = W*L
gap = 16u
f = 0.5 * eps * area * V^2 / gap^2 
f3d * [d][ F = f   oz=-pi/2] %actuator force

h = 200u
R1 = 1300u
R2 = 2000u
pi = 3.14159265358979
v = h * pi * (R1^2 - R2^2)
d = 2300
m = v*d
g = 9.8
f = m*g
%f3d * [b][ F=f/4 oy=pi/2] %wieght of the mass

