%net=cho_load('circringpos1.m');[q,k]=cho_dc(net);figure(1);cho_display(net);
%net=cho_load('circringpos1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'b','z'))
%net=cho_load('circringpos1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'b','z'));q50=q(1:6,1)
%[Q,khat,fhat]=matrixcondensation_endnode(k,[0;0;-200e-6;0;0;0],0);Q50=Q
   
uses mumps.net
uses stdlib.net
res=50
finalangle=pi/4
%finalangle=pi/2
%finalangle=pi
subnet arc[n(0) n(res)] [radius=100u]
[
  for k=0:res-1
  [
    pos    *      [o n(k)] [y=radius-radius*cos(k*finalangle/res) x=sin(k*finalangle/res)*radius]
    pos    *      [o n(k+1)] [y=radius-cos((k+1)*finalangle/res)*radius x=sin((k+1)*finalangle/res)*radius]
    beam3d parent [n(k) n(k+1)] [w=10u]
  ]
]

%anchor p1 [a] [l=5u w=5u]
%f3d * [b] [oy=pi/2 F=200u]
arc p1 [a b] []


