%net=cho_load('circringpos1.m');[q,k]=cho_dc(net);figure(1);cho_display(net);
%net=cho_load('circringpos1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'b','z'))
%net=cho_load('circringpos1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'b','z'));q50=q(1:6,1)
%[Q,khat,fhat]=matrixcondensation_endnode(k,[0;0;-200e-6;0;0;0],0);Q50=Q
   
uses mumps.net
uses stdlib.net
resolution=30

subnet semicircularbeam[n(0) n(resolution)] [radius=* finalangle=* w=* h=*]
[
  for k=0:resolution-1
  [
    pos    *      [o n(k)] [x=cos(k*finalangle/resolution)*radius y=sin(k*finalangle/resolution)*radius]
    pos    *      [o n(k+1)] [x=cos((k+1)*finalangle/resolution)*radius y=sin((k+1)*finalangle/resolution)*radius]
    beam3d parent [n(k) n(k+1)] [w=w h=h]
  ]
]

