
%net=cho_load('mirror9.m');[q,k]=cho_dc(net);figure(2);cho_display(net,q);shownodes(net);q(lookup_coord(net,'M','y')),q(lookup_coord(net,'M3','y'))
%net=cho_load('mirror9.m');figure(1);cho_display(net);
%net=cho_load('mirror9.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%clear all;net=cho_load('mirror3.m');figure(1);cho_display(net);shownodes(net);

uses mumps.net
uses comb.m
uses perferated.m
%uses predisplacedbeam2.m

%anchor p1 [a1][l=10u h=10u w=20u oz=pi]
%predisplacedbeam4 p1 [a1 a11][l=100u w=10u h=2u qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=55u qz2=0 qox2=0 qoy2=0 qoz2=0] 

%comb p1 [left middle right][V=500 numberfingers=11 gap=6u fingerwidth=2u fingerheight=2u fingerlength=20u supportwidth=2u supportheight=2u oz=pi/4] 
%anchor p1 [left][l=10u h=10u w=20u oz=pi]


wplate=700u-100u
Rplate=350u-50u
hplate=30u
arcplate=pi/2
ozplate=pi/2
semicircularbeam p1 [p0 p1] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=0*ozplate]
beam3d p1 [p p0][l=100u w=20u h=20u oz=0]
anchor p1 [p][l=40u h=40u w=40u oz=pi]
beam3dlink p1 [p1 b][l=200u w=20u h=20u oz=0 L1=1000u]

semicircularbeam p1 [p1 b] [radius=200u w=400u h=hplate alpha=arcplate L1=1000u oz=pi/2]

f3d * [b][F=1000u oy=pi/2]

%semicircularbeam p1 [p0 p1] [radius=Rrim w=wrim h=hrim alpha=arcplate oz=0*ozplate oy1=pi/2 oz1=(-pi/2+theta) L1=wplate/2+wrim/2 oy2=pi oz2=theta L2=wplate/2+wrim/2]



%beam3d p1 [p p0][l=100u w=10u h=2u oz=0 ]
%beam3d p1 [p0 m][l=100u w=10u h=10u oz=0 ]
%beam3d p1 [m p1][l=100u w=10u h=10u oz=pi/2 ]
%beam3dlink p1 [p0 O][l=200u w=2u h=2u oz=0 L1=100u oz1=-pi/2]
%beam3d p1 [p0 M2][l=100u w=20u h=2u oz=-pi/2 ]
%beam3d p1 [M2 M3][l=200u w=20u h=2u oz=0 ]
%beam3dlink p1 [O p1][l=200u w=2u h=2u L2=100u oz2=pi/2 oz=pi/2]
%anchor p1 [p][l=40u h=4u w=40u oz=pi]
%f3d * [p1][F=2000u oy=0*pi/2]

%beam3d p1 [p1 X2][l=100u w=2u h=2u oz=0]
%beam3d p1 [X2 X3][l=200u w=2u h=2u oz=-pi/2 ]



%semicircularbeam p1 [p1 p2] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=1*ozplate]
%semicircularbeam p1 [p2 p3] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=2*ozplate]
%semicircularbeam p1 [p3 p4] [radius=Rplate w=wplate h=hplate alpha=arcplate oz=3*ozplate]
%rim
wrim=50u
Rrim=Rplate*2+25u
hrim=50u
theta=atan(((hrim-hplate)/2)/(wplate/2+wrim/2))
%semicircularbeam p1 [p0 p1] [radius=Rrim w=wrim h=hrim alpha=arcplate oz=0*ozplate oy1=pi/2 oz1=(-pi/2+theta) L1=wplate/2+wrim/2 oy2=pi oz2=theta L2=wplate/2+wrim/2]
%semicircularbeam p1 [p1 up2] [radius=Rrim w=wrim h=hrim alpha=arcplate oz=1*ozplate ]%oy1=pi/2 oz1=(-pi/2+theta) L1=wplate/2+wrim/2 oy2=pi oz2=theta L2=wplate/2+wrim/2]
%semicircularbeam p1 [p2 up3] [radius=Rrim w=wrim h=hrim alpha=arcplate oz=2*ozplate ]%oy1=pi/2 oz1=(-pi/2+theta) L1=wplate/2+wrim/2 oy2=pi oz2=theta L2=wplate/2+wrim/2]
%semicircularbeam p1 [p3 up4] [radius=Rrim w=wrim h=hrim alpha=arcplate oz=3*ozplate ]%oy1=pi/2 oz1=(-pi/2+theta) L1=wplate/2+wrim/2 oy2=pi oz2=theta L2=wplate/2+wrim/2]

%beam3d p1 [p1p3 a][l=wplate*2+200u w=20u oz=0]
%f3d * [p1][F=1000u oy=0*pi/2]

%anchor p1 [a][l=10u h=10u w=20u oz=pi]
%beam3d p1 [a p0][l=150u w=20u oz=pi/2]
%beam3d p1 [p3 e][l=150u w=20u oz=pi/2]
%f3d * [e][F=1u oy=pi/2]

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

