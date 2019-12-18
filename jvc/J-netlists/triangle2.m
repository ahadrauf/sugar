%net=cho_load('triangle2.m');figure(1);cho_display(net);

uses mumps2.net

res=10
anchor p1 [A][l=10u w=10u oz=-90]
theta=30
alpha=30
W=10u
L1=150u
L2=100u
beam3d p1 [A b(1,3)][w=W l=L1/res oz=(theta+alpha)] 
beam3d p1 [A a(1,3)][w=W l=L2/res oz=theta] 
pi=3.141592

for k=1:res [
   beam3d p1 [b(k,3) b(k+1,3)][w=W l=L1/res oz=(theta+alpha)] 
   beam3d p1 [a(k,3) a(k+1,3)][w=W l=L2/res oz=theta] 
   x1=(k+1)*L1/res*cos(theta*pi/180)
   y1=(k+1)*L1/res*sin(theta*pi/180)
   x2=(k+1)*L2/res*cos((theta+alpha)*pi/180)
   y2=(k+1)*L2/res*sin((theta+alpha)*pi/180)
   L23=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
   phi=(pi/2-atan((y1-y2)/(x1-x2)))*180/pi
   W12=L2/res 
   beam3d p1 [a(k+1,3) b(k+1,3)][w=W12 l=L23 oz=phi] 
   for j=1:res [
         r=j
         ]
] 

f3d * [a(10,3)][M=1u oy=-90] 


