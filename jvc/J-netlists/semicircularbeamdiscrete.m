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

%discreteness
resolution=20 %number of straight beams that curved beam is broken into.
%geometry
param alpha=pi/4, radius=100u, w=10u, h=2u, ox=0, oy=0, oz=0
len = sqrt((radius*sin(alpha/resolution))^2 + (radius - radius*cos(alpha/resolution))^2 )

subnet semicircularbeamdiscrete [n(0) n(resolution)] []
[
  for k=0:resolution-1
  [
    pos    *      [o n(k)]   [y=radius-radius*cos(k*alpha/resolution)     x=sin(k*alpha/resolution)*radius]
    pos    *      [o n(k+1)] [y=radius-cos((k+1)*alpha/resolution)*radius x=sin((k+1)*alpha/resolution)*radius]
    beam3d parent [n(k) n(k+1)] [l=len w=w h=h]
  ]
]

semicircularbeamdiscrete layer1 [node1 node2] [ox=ox oy=oy oz=oz]

