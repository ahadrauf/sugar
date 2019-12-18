%micromirror
%net=cho_load('mirror10.m');figure(1);cho_display(net);
%net=cho_load('mirror10.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%clear all;net=cho_load('mirror3.m');figure(1);cho_display(net);shownodes(net);
%clear all;tic;net=cho_load('mirror10.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);toc

% L I B R A R Y  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uses mumps.net
uses comb.m
uses perferated.m
%uses scs100.net   <--- make this
%uses predisplacedbeam2.m   <--- check node positioning

lowSCS=10u
highSCS=50u
upperSCS=15u

% P L A T E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plate parameters
w_plate=300u
radius_plate=w_plate/2
h_plate=upperSCS
arcplate=pi/2
oz_plate=pi/2

%rim parameters
w_rim=50u
radius_rim=w_plate+w_rim/2
h_rim=highSCS
theta=atan((h_rim/2-h_plate/2)/(w_plate/2+w_rim/2))
L=sqrt((w_plate/2+w_rim/2)^2+(h_rim/2-h_plate/2)^2)
L1=L
L2=L

%plate
semicircularbeam p1 [plate_bottom plate_right]  [radius=radius_plate w=w_plate h=h_plate alpha=arcplate oz=0*oz_plate]
semicircularbeam p1 [plate_right  plate_top]    [radius=radius_plate w=w_plate h=h_plate alpha=arcplate oz=1*oz_plate]
semicircularbeam p1 [plate_top    plate_left]   [radius=radius_plate w=w_plate h=h_plate alpha=arcplate oz=2*oz_plate]
semicircularbeam p1 [plate_left   plate_bottom] [radius=radius_plate w=w_plate h=h_plate alpha=arcplate oz=3*oz_plate]

%rim
semicircularbeam p1 [plate_bottom plate_right]  [radius=radius_rim w=w_rim h=h_rim alpha=arcplate oz=0*oz_plate oy1=pi/2 oz1=(-pi/2+theta) L1=L oy2=theta oz2=pi L2=L] 
semicircularbeam p1 [plate_right  plate_top]    [radius=radius_rim w=w_rim h=h_rim alpha=arcplate oz=1*oz_plate oy1=pi/2 oz1=(-pi/2+theta) L1=L oy2=theta oz2=pi L2=L]
semicircularbeam p1 [plate_top    plate_left]   [radius=radius_rim w=w_rim h=h_rim alpha=arcplate oz=2*oz_plate oy1=pi/2 oz1=(-pi/2+theta) L1=L oy2=theta oz2=pi L2=L]
semicircularbeam p1 [plate_left   plate_bottom] [radius=radius_rim w=w_rim h=h_rim alpha=arcplate oz=3*oz_plate oy1=pi/2 oz1=(-pi/2+theta) L1=L oy2=theta oz2=pi L2=L]

% H I N G E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hinge parameters
w_hinge2plate=20u
w_hinge=9u
w_hinge=9u
h_hinge=highSCS
h_hinge2plate=highSCS
l_hingearm=340u-w_hinge2plate/2-w_hinge/2
l_hingecenter=30u/2
l_hinge2plate=60u

%tether parameters
l_tether=735u
h_tether=lowSCS
w_tether=12u

%beam from hinge to plate
beam3dlinkcorner p1 [hinge_right plate_left][l=l_hinge2plate w=w_hinge2plate h=h_hinge2plate L1=w_hinge/2 L2=w_plate/2+w_rim] %L2 is rigid link to plate

%hinge
beam3dlink p1 [hinge_right  hinge_rightcorner] [l=l_hingearm    w=w_hinge h=h_hinge oz=pi/2 L1=w_hinge2plate/2 L2=w_hinge/2] %Right arm
beam3dlink p1 [hinge_center hinge_rightcorner] [l=l_hingecenter w=w_hinge h=h_hinge L2=w_hinge/2] %Center right
beam3dlink p1 [hinge_leftcorner  hinge_center] [l=l_hingecenter w=w_hinge h=h_hinge L1=w_hinge/2] %Center left
beam3dlink p1 [hinge_leftcorner    hinge_left] [l=l_hingearm    w=w_hinge h=h_hinge oz=-pi/2 L1=w_hinge/2 L2=w_hinge2plate/2] %Left arm

%tether for the hinge
beam3dlink p1 [hinge_center bridgeleft][l=l_tether w=w_tether h=h_tether oz=-pi/2 L1=h_hinge/2-h_tether/2 oy1=pi/2] 

% P E R F O R A T E D   A R M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_perf2hinge=40u
w_perf2hinge=20u
h_perf2hinge=upperSCS

%beam from perforated to hinge
beam3dlinkconrner p1 [perf_right hinge_left][l=l_perf2hinge w=w_perf2hinge h=h_perf2hinge L2=w_hinge/2]

%perforated beam
PerferatedW=20u
PerferatedH=15u
PerferatedL=400u
PerferatedHoles=24
HoleVertW=20u-12u*2
HoleHoriW=4.6u
perferated p1 [perf_left perf_right][numberholes=PerferatedHoles w=PerferatedW h=PerferatedH l=PerferatedL whorizontal=HoleHoriW wvertical=HoleVertW] 

%anchor for perforated beam
anchor p1 [perf_left][l=40u h=40u w=40u oz=pi] 

% B R I D G E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bridgeL=860u/2
bridgeH=50u
bridgeW=30u
beam3dlink p1 [bridgeleft bridgecenter][l=bridgeL w=bridgeW h=bridgeH] 






%f3d * [plate_right][F=100u oy=-pi/2]

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

%anchor p1 [a1][l=10u h=10u w=20u oz=pi] 
%predisplacedbeam4 p1 [a1 a11][l=100u w=10u h=2u qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=55u qz2=0 qox2=0 qoy2=0 qoz2=0] 
%comb p1 [left middle right][V=500 numberfingers=11 gap=6u fingerwidth=2u fingerheight=2u fingerlength=20u supportwidth=2u supportheight=2u oz=pi/4]  
%anchor p1 [left][l=10u h=10u w=20u oz=pi] 

