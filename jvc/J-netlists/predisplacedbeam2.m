uses mumps.net
resolution=20 %number of straight beams that curved beam is broken into.
param alpha=pi/4, radius=100u, w=10u, h=2u, ox=0, oy=0, oz=0

subnet semicircularbeamdiscrete [n(0) n(resolution)] []
[
  for k=0:resolution-1
  [
    pos    *      [o n(k)]   [y=radius-radius*cos(k*alpha/resolution)     x=sin(k*alpha/resolution)*radius]
    pos    *      [o n(k+1)] [y=radius-cos((k+1)*alpha/resolution)*radius x=sin((k+1)*alpha/resolution)*radius]
    beam3d parent [n(k) n(k+1)] [w=w h=h l=10u]
  ]
]

semicircularbeamdiscrete p1 [node1 node2] [ox=ox oy=oy oz=oz]

