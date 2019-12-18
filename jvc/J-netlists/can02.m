% net=cho_load('can00.m');[q]=cho_dc(net);

uses mumpsx.net
param p
anchor  p1 [A]   [l=10u w=10u]
beam3d  p3 [A B] [l=200u w=4u h=2u oz=0 t=p*p']
f3d     *  [B]   [F=1 oz=90]

