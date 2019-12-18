%net=cho_load('mirror7.m');figure(1);cho_display(net);
%net=cho_load('mirror7.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%clear all;net=cho_load('mirror3.m');figure(1);cho_display(net);shownodes(net);

uses mumps.net
uses comb.m
uses perferated.m
%uses predisplacedbeam2.m

%anchor p1 [a1][l=10u h=10u w=20u oz=pi]
%predisplacedbeam4 p1 [a1 a11][l=100u w=10u h=2u qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=55u qz2=0 qox2=0 qoy2=0 qoz2=0] 

%comb p1 [left middle right][V=500 numberfingers=11 gap=6u fingerwidth=2u fingerheight=2u fingerlength=20u supportwidth=2u supportheight=2u oz=pi/4] 
%anchor p1 [left][l=10u h=10u w=20u oz=pi]


wplate=700u
Rplate=350u
hplate=30u
arcplate=pi/4
ozplate=pi/4
anchor p1 [a][l=10u h=10u w=20u oz=pi]
beam3d p1 [a p0][l=150u w=20u oz=pi/2]
beam3d p1 [p3 e][l=150u w=20u oz=pi/2]
f3d * [e][F=1u oy=pi/2]
semicircularbeam p1 [p0 p1] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=0*ozplate L1=wplate/2 oz1=pi/2 ]
semicircularbeam p1 [p1 p2] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=1*ozplate ]
semicircularbeam p1 [p2 p3] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=2*ozplate L2=wplate/2 oz2=-pi/4]

%circular plate
wplate=700u
Rplate=350u
hplate=30u
arcplate=pi/4
ozplate=pi/4
%semicircularbeam p1 [p0 p1] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=0*ozplate L1=wplate/2 oz1=pi/2 L2=wplate/2 oz2=-pi/4]
%semicircularbeam p1 [p1 p2] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=1*ozplate L1=wplate/2 oz1=pi/2 L2=wplate/2 oz2=-pi/4]
%semicircularbeam p1 [p2 p3] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=2*ozplate L1=wplate/2 oz1=pi/2 L2=wplate/2 oz2=-pi/4]
%semicircularbeam p1 [p3 p4] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=3*ozplate L1=wplate/2 oz1=pi/2 L2=wplate/2 oz2=-pi/4]
%semicircularbeam p1 [p4 p5] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=4*ozplate L1=wplate/2 oz1=pi/2 L2=wplate/2 oz2=-pi/4]
%semicircularbeam p1 [p5 p6] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=5*ozplate L1=wplate/2 oz1=pi/2 L2=wplate/2 oz2=-pi/4]
%semicircularbeam p1 [p6 p7] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=6*ozplate L1=wplate/2 oz1=pi/2 L2=wplate/2 oz2=-pi/4]
%semicircularbeam p1 [p7 p0] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=7*ozplate L1=wplate/2 oz1=pi/2 L2=wplate/2 oz2=-pi/4]

%anchor p1 [a][l=10u h=10u w=20u oz=pi]
%beam3dlinkcorner p1 [a b][l=100u w=10u h=2u L2=50u oz=pi/4] 
%beam3dlinkcorner p1 [b c][l=100u w=10u h=2u L2=0 oz=pi/2] 
%beam3dlinkcorner p1 [c d][l=100u w=10u h=2u L1=0u] 
%f3d * [c][F=1000u oz=0]

%beam3dlink p1 [a b][l=100u w=20u h=2u oz=-pi/2 L1=supportwidth/2] 
%beam3d p1 [b4 B][l=150u w=20u oz=pi/2]

%f3d * [B][F=10u oy=pi/2]

%anchor p1 [left][l=10u w=10u h=10u oz=-pi]
%perferated p1 [left right][numberholes=4 w=20u h=10u l=50u whorizontal=2u wvertical=2u] 

