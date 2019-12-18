%predisplacedbeam2.m
%This netlist makes a predisplaced beam that starts off growing along the x-axis and deflects due to user defined displacements
%at node1 and node2. Only angle displacements are used on node1. Both translational and angular displacements on node2.
%Useful for generating the 12*12 effective stiffness matrix for a predisplaced beam element.
%Running commands:
%net=cho_load('predisplaced2.m');[q,k]=cho_dc(net);figure(1);cho_display(net);
%By: Jason Vaughn Clark, Nov2001

%libraries   
uses mumps.net
%parameterization
resolution=10
param l=100u, w=10u, h=2u, qx1=0, qy1=0, qox1=0, qoy1=0, qoz1=0, qx2=0, qy2=0, qz2=0, qox2=0, qoy2=0, qoz2=0 
param ox=0, oy=0, oz=0

%subnet predisplacedbeam4 [a(1) a(resolution+1)][l=* w=* h=* qx1=* qy1=* qox1=* qoy1=* qoz1=* qx2=* qy2=* qz2=* qox2=* qoy2=* qoz2=* ]

subnet predisplacedbeam4 [a(1) a(10+1)][]
[
L=l %Original beam length.
%These are not needed for now.
qx1=0
qy1=0
qz1=0
%Hermite polynomial coefficients
yf0=qy1
yf00=qoz1
yfL=qy2
yfLL=qoz2
yf0L=(yfL-yf0)/L
yf00L=(yf0L-yf00)/L
yf0LL=(yfLL-yf0L)/L
yf00LL=(yf0LL-yf00L)/L
zf0=qz1
zf00=qoy1
zfL=-qz2%
zfLL=-qoy2
zf0L=(zfL-zf0)/L
zf00L=(zf0L-zf00)/L
zf0LL=(zfLL-zf0L)/L
zf00LL=(zf0LL-zf00L)/L
%discretized beams
for i=1:resolution
[
   s1=L*(i-1)/resolution %first sample point
   s2=L*i/resolution %second sample point
   hx1=qx1+(1+(qx2-qx1)/L)*s1 
   hy1=yf0+s1*(yf00+s1*(yf00L+(s1-L)*yf00LL))
   hz1=zf0+s1*(zf00+s1*(zf00L+(s1-L)*zf00LL))
   hx2=qx1+(1+(qx2-qx1)/L)*s2 
   hy2=yf0+s2*(yf00+s2*(yf00L+(s2-L)*yf00LL))
   hz2=zf0+s2*(zf00+s2*(zf00L+(s2-L)*zf00LL))
   %slope from sample point 1 to 2   
   oy=atan((hz2-hz1)/(hx2-hx1))
   oz=atan((hy2-hy1)/(hx2-hx1))
   ox=qox1*(resolution-i)/(resolution-1)+qox2*(i-1)/(resolution-1)
   %discretized beam   
   beam3d parent [a(i) a(i+1)][l=(L+qx2-qx1)/resolution w=w h=h ox=ox oy=oy oz=oz]    
]%for i
]%subnet

predisplacedbeam4 p1 [node1 node2][ox=ox oy=oy oz=oz]

%predisplacedbeam4 p1 [node1 node2][l=l w=w h=h qx1=qx1 qy1=qy1 qox1=qox1 qoy1=qoy1 qoz1=qoz1 qx2=qx2 qy2=qy2 qz2=qz2 qox2=qox2 qoy2=qoy2 qoz2=qoz2]



