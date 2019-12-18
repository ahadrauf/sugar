%micromirror
%net=cho_load('mirror13.m');figure(1);cho_display(net);
%net=cho_load('mirror13.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%clear all;tic;net=cho_load('mirror13.m');figure(1);cho_display(net);toc
%clear all;tic;net=cho_load('mirror13.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);toc

% L I B R A R Y  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uses mumps.net
uses comb.m
uses perforated.m
uses perforatedcomb.m
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
h_hinge=highSCS
h_hinge2plate=highSCS
l_hingearm=340u-w_hinge2plate/2-w_hinge/2
l_hingecenter=30u/2
l_hinge2plate=60u

%tether parameters
l_tether=735u
h_tether=lowSCS
w_tether=12u

theta_hinge2plate=atan((h_hinge2plate/2-h_plate/2)/(w_plate/2+w_rim))
L2_hinge2plate=sqrt((w_plate/2+w_rim)^2+(h_hinge2plate/2)^2)

subnet hinge [plate_left hinge_center hinge_left][]
[
%beam from hinge to plate
beam3dlink p1 [hinge_right plate_left][l=l_hinge2plate w=w_hinge2plate h=h_hinge2plate L1=w_hinge/2 L2=L2_hinge2plate oy2=-theta_hinge2plate] %L2 is rigid link to plate

%hinge
beam3dlinkcorner p1 [hinge_right  hinge_rightcorner] [l=l_hingearm    w=w_hinge h=h_hinge oz=pi/2 L1=w_hinge2plate/2 L2=w_hinge/2] %Right arm
beam3dlink p1 [hinge_center hinge_rightcorner] [l=l_hingecenter w=w_hinge h=h_hinge L2=w_hinge/2] %Center right
beam3dlink p1 [hinge_leftcorner  hinge_center] [l=l_hingecenter w=w_hinge h=h_hinge L1=w_hinge/2] %Center left
beam3dlinkcorner p1 [hinge_leftcorner    hinge_left] [l=l_hingearm    w=w_hinge h=h_hinge oz=-pi/2 L1=w_hinge/2 L2=w_hinge2plate/2] %Left arm
]

hinge p1 [plate_left hinge_center hinge_left][]

%tether for the hinge
beam3dlink p1 [hinge_center bridgeleft][l=l_tether w=w_tether h=h_tether oz=-pi/2 L1=h_hinge/2-h_tether/2 oy1=pi/2 L2=h_hinge/2-h_tether/2 oy2=-pi/2] 

% P E R F O R A T E D   A R M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%perforated beam
l_perf1=71u
l_perf2=400u-l_perf1
w_perf=20u
w_holevert=4.6u
h_perf=upperSCS
perfholes=24
w_holehori=(w_perf-11u)/2

l_perf2hinge=40u
w_perf2hinge=20u
h_perf2hinge=highSCS

h_perfoverhang=highSCS-upperSCS
w_perfoverhang=w_holehori
l_perfoverhang=l_perf1
L_perfoverhang=sqrt((w_perf/2-w_perfoverhang/2)^2+(h_perf/2+h_perfoverhang/2)^2)
theta_perfoverhang=atan((w_perf/2-w_perfoverhang/2)/(h_perf/2+h_perfoverhang/2))

subnet perferatedbeam [hinge_left][]
[
%beam from perforated to hinge
beam3dlink p1 [perf_right hinge_left][l=l_perf2hinge w=w_perf2hinge h=h_perf2hinge L2=w_hinge/2 L1=(h_perf2hinge/2-h_perf/2) oy1=pi/2]

%perforated beam
perforated p1 [perf_middle perf_right][numberholes=perfholes-20 w=w_perf h=h_perf l=l_perf1 whorizontal=w_holehori wvertical=w_holevert] 
perforated p1 [perf_left  perf_middle][numberholes=perfholes-4  w=w_perf h=h_perf l=l_perf2 whorizontal=w_holehori wvertical=w_holevert] 

%overhang of perferated
beam3dlink p1 [perf_middle perf_right][l=l_perfoverhang h=h_perfoverhang w=w_perfoverhang L1=L_perfoverhang oz1=-theta_perfoverhang oy1=pi/2 L2=L_perfoverhang oz2=theta_perfoverhang oy2=-pi/2] 
beam3dlink p1 [perf_middle perf_right][l=l_perfoverhang h=h_perfoverhang w=w_perfoverhang L1=L_perfoverhang oz1=theta_perfoverhang oy1=pi/2 L2=L_perfoverhang oz2=-theta_perfoverhang oy2=-pi/2] 

%anchor for perforated beam
anchor p1 [perf_left][l=40u h=40u w=40u oz=pi]  
]

perferatedbeam p1 [hinge_left][]

% B R I D G E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_bridge=l_hingecenter+w_hinge+l_hinge2plate+w_rim+w_plate+w_tether/2
h_bridge=50u
w_bridge=30u
beam3dlink p1 [bridgeleft bridgecenter][l=l_bridge w=w_bridge h=h_bridge L1=-w_tether/2] 

% C O M B D R I V E   A R R A Y  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

perfholes=6
w_perf=30u
h_perf=highSCS
w_horizontal=10u
w_vertical=8u
l_perf=100u
perforated p1 [bridgecenter comb1][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 

subnet combarray [comb1 comb4][]
[
w_finger=3.2u
l_finger=30u
gap=3u
numberfingers=70
w_perf=30u
whorizontal=(w_perf-12u)/2
h_perf=highSCS
w_vertical=6u

perforated p1 [bridgecenter comb1][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 
perforatedcomb p1 [leftcomb1  comb1][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical] 
perforatedcomb p1 [comb1 rightcomb1][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical] 

l_perf=140u
perforated p1     [comb1      comb2][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 
perforatedcomb p1 [leftcomb2  comb2][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical] 
perforatedcomb p1 [comb2 rightcomb2][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical] 

perforated p1     [comb2      comb3][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 
perforatedcomb p1 [leftcomb3  comb3][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical] 
perforatedcomb p1 [comb3 rightcomb3][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical] 

perforated p1     [comb3      comb4][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 
perforatedcomb p1 [leftcomb4  comb4][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical] 
perforatedcomb p1 [comb4 rightcomb4][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical] 
]

combarray p1 [comb1 comb4][]

l_perf=140u
perforated p1 [comb4 comb5][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 



% P R E D I S P L A C E D   B E A M S  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_predisp=1100u
w_predisp=6.5u
h_predisp=highSCS
y2_disp=-25u
displink=30u
predisplacedbeam4 p1 [bridgecenter disp1][l=l_predisp w=w_predisp h=h_predisp qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=y2_disp qz2=0 qox2=0 qoy2=0 qoz2=0 L1=displink oz1=-pi/2] 
predisplacedbeam4 p1 [comb5        disp2][l=l_predisp w=w_predisp h=h_predisp qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=y2_disp qz2=0 qox2=0 qoy2=0 qoz2=0] 

anchor p1 [disp1][l=40u h=40u w=40u] 

%f3d * [B][F=10u oy=pi/2]

%anchor p1 [left][l=10u w=10u h=10u oz=-pi]
%perforated p1 [left right][numberholes=4 w=20u h=10u l=50u whorizontal=2u wvertical=2u] 

%anchor p1 [a1][l=10u h=10u w=20u oz=pi] 
%predisplacedbeam4 p1 [a1 a11][l=100u w=10u h=2u qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=55u qz2=0 qox2=0 qoy2=0 qoz2=0] 
%comb p1 [left middle right][V=500 numberfingers=11 gap=6u fingerwidth=2u fingerheight=2u fingerlength=20u supportwidth=2u supportheight=2u oz=pi/4]  
%anchor p1 [left][l=10u h=10u w=20u oz=pi] 

