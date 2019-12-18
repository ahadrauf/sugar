%semicircularbeam.m
%This netlist makes a semicircular beam that starts off growing along the x-axis, curving up in quadrant-1.
%Useful for generating the 12*12 effective stiffness matrix for a semicircular beam element.
%The parameters are alpha (total angle swept), radius (radius from origin to neutral beam axis), w (beam width), h (layer thickness)
%   ox oy oz (orientations about the x y and z-axes).
%Running commands:
%net=cho_load('semicircularbeamdiscrete.m');[q,k]=cho_dc(net);figure(1);cho_display(net);
%By: Jason Vaughn Clark, Nov2001

%libraries   
uses generalprocess.m

%geometry
param l=100u, w=10u, h=2u, qx1=0, qy1=0, qox1=0, qoy1=0, qoz1=0, qx2=0, qy2=0, qz2=0, qox2=0, qoy2=0, qoz2=0
param ox=0, oy=0, oz=0
param resolution=50 

subnet predisplacedbeamdiscrete [n(0) n(resolution)] []
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

   for k=0:resolution-1
   [
      s1=L*(k-1)/resolution %first sample point
      s2=L*k/resolution %second sample point
      hx1=qx1+(1+(qx2-qx1)/L)*s1 
      hy1=yf0+s1*(yf00+s1*(yf00L+(s1-L)*yf00LL))
      hz1=zf0+s1*(zf00+s1*(zf00L+(s1-L)*zf00LL))
      hx2=qx1+(1+(qx2-qx1)/L)*s2 
      hy2=yf0+s2*(yf00+s2*(yf00L+(s2-L)*yf00LL))
      hz2=zf0+s2*(zf00+s2*(zf00L+(s2-L)*zf00LL))
      %slope from sample point 1 to 2   
      oy=atan((hz2-hz1)/(hx2-hx1))
      oz=atan((hy2-hy1)/(hx2-hx1))
      ox=qox1*(resolution-k)/(resolution-1)+qox2*(k-1)/(resolution-1)
      %discretized beam   
      beam3d parent [n(k) n(k+1)][l=(L+0*(qx2-qx1))/resolution w=w h=h ox=ox oy=oy oz=oz resolution=resolution qx1=qx1 qy1=qy1 qox1=qox1 qoy1=qoy1 qoz1=qoz1 qx2=qx2 qy2=qy2 qz2=qz2 qox2=qox2 qoy2=qoy2 qoz2=qoz2]    
   ]
]

predisplacedbeamdiscrete layer1 [node1 node2] [ox=ox oy=oy oz=oz w=w h=h ox=ox oy=oy oz=oz resolution=resolution qx1=qx1 qy1=qy1 qox1=qox1 qoy1=qoy1 qoz1=qoz1 qx2=qx2 qy2=qy2 qz2=qz2 qox2=qox2 qoy2=qoy2 qoz2=qoz2]

