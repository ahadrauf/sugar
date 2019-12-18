% net=cho_load('semicirc1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q); q(2)

uses mumps.net

res = 12
theta = pi/2/res
R = 97.5u
L = sqrt( (R*sin(theta))^2 + (R - R*cos(theta))^2 )
W = 5u
H = 2u


anchor   p1 [a(0)]  [l=10u w=10u h=10u oz=-pi/2]

n=1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(0) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 
n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 
n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 

n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 
n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 

n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 


n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 
n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 
n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 

n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 
n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 
n=n+1
d = R*cos((n-1)*theta) - R*cos(n*theta)
phi = asin(d/L)
oz = pi/2 + phi
beam3d p1 [a(n-1) a(n)][w=W h=H l=L Youngsmodulus=170e9 oz=oz ] 

f3d * [a(n)]  [F=7u oz=0 ox=0 oy=-pi/2] 

