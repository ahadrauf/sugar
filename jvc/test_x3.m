
res=10
theta=30
alpha=30
W=10
L1=100
L2=100
%beam3d p1 [A b(1)][w=W l=L1/res oz=(theta+alpha)] 
%beam3d p1 [A a(1)][w=W l=L2/res oz=theta] 

for k=1:2
%   beam3d p1 [b(k) b(k+1)][w=W l=L1/res oz=(theta+alpha)] 
%   beam3d p1 [a(k) a(k+1)][w=W l=L2/res oz=theta] 
   x1=k*L1/res*cos(theta)
   y1=k*L1/res*sin(theta)
   x2=k*L2/res*cos(theta+alpha)
   y2=k*L2/res*sin(theta+alpha)
   L23=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
   phi=90-180/3.1456*atan((y2-y1)/(x2-x1))
   W12=L2/res 
%   beam3d p1 [a(k+1) b(k+1)][w=W12 l=L23  oz=phi] 
end
