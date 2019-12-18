

L = length_guess
K = desired_stiffness
for w=w1:setp:w2
   
   p.f=f;net=cho_load('nlvslin.m',p);[q]=cho_dc(net);figure(1);cho_display(net,q);q(2),end

uses mumpsx.net
anchor    p1 [A]   [l=10u w=10u] 

beam3dnl5 p1 [A B] [l=200u w=2u h=2u oz=0] %Nonlinear beam
beam3d    p1 [B b] [l=4u w=2u h=2u oz=90] %beam tip

beam3d    p1 [A C] [l=200u w=2u h=2u oz=0] %Linear beam
beam3d    p1 [C c] [l=4u w=2u h=2u oz=90] %beam tip

param f
f3d        * [B]   [F=f oz=90] 
f3d        * [C]   [F=f oz=90] 

