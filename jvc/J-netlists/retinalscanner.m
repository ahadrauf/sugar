%net=cho_load('retinalscanner.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);

%uses generalprocess.m
%uses comb.m
%uses perforated.m
%uses perforatedcomb4.m

uses mumps.net
h=4u

w_anchor=200u
l_anchor=300u
w_anchor_attach=13u
l_anchor_attach=35u

w_ubeam=9u
l_ubeam=300u
l_ubeam2=l_ubeam + 30u
r_ubeam=22u



%anchor 
anchor      p1 [b]  [l=1u w=1u h=1u] %anchor
%beam3d      p1 [a1 a2]  [l=l_anchor w=w_anchor h=h ] %anchor
%beam3d  p1 [a2 b]   [l=l_anchor_attach w=w_anchor_attach h=h       L1=w_anchor/2-w_anchor_attach/2 oz1=-pi/2 L2=-w_ubeam/2] %anchor attach

%U-shaped spring 
beam3d       p1 [b c]    [l=l_ubeam w=w_ubeam h=h oz=pi/2       L1=w_anchor_attach/2 ] %left vertical
%semicircularbeam p1 [c d]    [w=w_ubeam h=h radius=r_ubeam alpha=pi/2 oz=pi/2 ox=pi] %left curve
%semicircularbeam p1 [d e]    [w=w_ubeam h=h radius=r_ubeam alpha=pi/2 oz=0 ox=pi] %left curve
%beam3dlink       p1 [d f1]    [l=l_ubeam2 w=w_ubeam h=h oz=-pi/2     L2=w_anchor_attach/2 ] %left vertical
%beam3dlink       p1 [f1 f2 ]    [l=l_ubeam2 w=w_ubeam h=h oz=-pi/2     L2=w_anchor_attach/2 ] %left vertical
%beam3dlink       p1 [f2 f3 ]    [l=l_ubeam2 w=w_ubeam h=h oz=-pi/2     L2=w_anchor_attach/2 ] %left vertical
%f3d * [f][F=100u]

%l_scanner_attach=61u
%w_scanner_attach=w_anchor_attach
%beam3dlink       p1 [g f]    [l=l_scanner_attach w=w_scanner_attach h=h  L2=-w_ubeam/2 ] 

%l_scanner=(300u / 2) - w_scanner_attach
%w_scanner=18u
%beam3dlink       p1 [gm g]    [l=l_scanner w=w_scanner_attach h=h oz=pi/2 L2=w_scanner_attach/2 ] 

beam3d p1 [c e]    [l=r_ubeam*2 w=w_ubeam h=h ] 
beam3d p1 [e e1]    [l=l_ubeam2 w=w_ubeam h=h oz=-pi/2] 
beam3d p1 [e1 e2]    [l=r_ubeam*2 w=w_ubeam h=h oz=pi] 
beam3d p1 [e2 e3]    [l=l_ubeam*2 w=w_ubeam/3 h=h oz=pi/2] 
f3d * [e1][F=100u]
