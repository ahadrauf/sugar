%running commands
%net=cho_load('semicircularbeamdiscrete.m');[q,k]=cho_dc(net);figure(1);cho_display(net);
%net=cho_load('circringpos1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'b','z'))
%net=cho_load('circringpos1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'b','z'));q50=q(1:6,1)
%[Q,khat,fhat]=matrixcondensation_endnode(k,[0;0;-200e-6;0;0;0],0);Q50=Q
   
%libraries   
uses mumps.net
uses stdlib.net

%discreteness
resolution=30 
%orientation
param ox=0
param oy=0
param oz=0
%geometry
param radius=100u
param finalangle=pi/4 %alpha
param w=10u
param h=2u 

subnet semicircularbeamdiscrete[n(0) n(resolution)] []
[
  for k=0:resolution-1
  [
    pos    *      [o n(k)] [y=radius-radius*cos(k*finalangle/resolution) x=sin(k*finalangle/resolution)*radius]
    pos    *      [o n(k+1)] [y=radius-cos((k+1)*finalangle/resolution)*radius x=sin((k+1)*finalangle/resolution)*radius]
    beam3d parent [n(k) n(k+1)] [w=w h=h]
  ]
]

semicircularbeamdiscrete p1 [a b][]


