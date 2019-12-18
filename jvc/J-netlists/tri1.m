% net=cho_load('tri1.m');[q]=cho_dc(net);cho_display(net,q);

uses mumpsx.net

anchor      p1 [A]      [l=10u w=10u h=10u oy=pi]
beam3d      p1 [A D]    [l=200u w=20u h=2u oz=45]
triangle2   p1 [A B C]  [l1=200u l2=200u angle=pi/4 h=2u oz=0]

%beam3d   p1 [B E]  [l=100u w=10u oz=0]
%f3d      *  [C]    [F=100000u oy=90]
%f3d      *  [C]    [M=1u oy=90]
%f3d      *  [E]    [F=3000u oz=90]
%f3d      *  [C]    [F=30000u oy=90]
