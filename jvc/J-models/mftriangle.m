%p.angle=90;p.l2=150e-6;net=cho_load('mftriangle.m',p);figure(1);cho_display(net);
%net=cho_load('mftriangle.m');figure(1);cho_display(net);
uses mumps2.net
%anchor p1 [A][l=10u w=10u oz=-90]
pi=3.141592
param ox=0, oy=0, oz=0
%param Poisson=0.3
%param thermcond
%param viscosity
%param fluid
%param density
%param Youngsmodulus
%param permittivity
%param sheetresistance
%param stress
%param straingradient
%param thermalexpansion
%param ambienttemperature
%param h
param l1=100u, l2=100u, angle=90, res=10
resolution=res
x2=l2*cos(angle*pi/180)
x1=l1
y2=l2*sin(angle*pi/180)
W=sqrt((x2+x1)*(x2+x1)+y2*y2)/2/(resolution)
oz=oz*180/pi
beam3d p1 [A b(1)][w=W l=l2/resolution oz=angle] 
beam3d p1 [A a(1)][w=W l=l1/resolution] 
%main diagonals
for k=1:resolution-1 [
   beam3d p1 [b(k) b(k+1)][w=W l=l2/resolution oz=angle] 
   beam3d p1 [a(k) a(k+1)][w=W l=l1/resolution] 
   x1=(k+1)*l1/resolution
   y1=0
   x2=(k+1)*l2/resolution*cos(angle*pi/180)
   y2=(k+1)*l2/resolution*sin(angle*pi/180)
   L21=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
   phi=-acos((x1-x2)/L21)*180/pi
   beam3d p1 [b(k+1) a(k+1)][w=W l=L21 oz=phi ] 
] 
%smallest diagonal
k=0
x1=(k+1)*l1/resolution
y1=0
x2=(k+1)*l2/resolution*cos(angle*pi/180)
y2=(k+1)*l2/resolution*sin(angle*pi/180)
L21=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
phi=-acos((x1-x2)/L21)*180/pi
beam3d p1 [b(k+1) a(k+1)][w=W l=L21 oz=phi t1=l1 t2=l2] 
%last diagonal
k=resolution
beam3d p1 [b(k) C][w=W l=l2/resolution oz=angle] 
beam3d p1 [a(k) B][w=W l=l1/resolution] 
x1=(k+1)*l1/resolution
y1=0
x2=(k+1)*l2/resolution*cos(angle*pi/180)
y2=(k+1)*l2/resolution*sin(angle*pi/180)
L21=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
phi=-acos((x1-x2)/L21)*180/pi
beam3d p1 [C B][w=W l=L21 oz=phi ] 
   
   