% theta=q(lookup_coord(net,'plate_top','rx')),cosine=q(lookup_coord(net,'bridgecenter','y'))

%verification:
%Mirror tilts about 0.157 rad for a cosine-beam displacement of 25um=2.5e-5m.

%micromirror
%net=cho_load('mirror20.m');figure(1);cho_display(net);
%net=cho_load('mirror20.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%clear all;tic;net=cho_load('mirror20.m');figure(1);cho_display(net);toc
%clear all;tic;net=cho_load('mirror20.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);toc

%f3d * [bridgecenter][F=1500u oz=-pi/2]

f3d * [bridgecenter][F=2300u oz=-pi/2]


% L I B R A R Y  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uses generalprocess.m
uses comb.m
uses perforated.m
uses perforatedcomb4.m
%uses predisplacedbeam2.m   <--- check node positioning

lowSCS=8u %8u:12u
highSCS=55u %50u
upperSCS=12u %15u

% P L A T E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plate parameters
w_plate=350u
radius_plate=w_plate/2
h_plate=upperSCS
arcplate=pi/2
oz_plate=pi/2

%rim parameters
w_rim=50u
radius_rim=w_plate-w_rim/2
h_rim=highSCS-upperSCS
theta_plate=atan((h_rim/2+h_plate/2)/(w_plate/2-w_rim/2))
L_plate=sqrt((w_plate/2-w_rim/2)^2+(h_rim/2+h_plate/2)^2)
L1_plate=L_plate
L2_plate=L_plate

subnet plate [plate_left plate_right plate_top][]
[
   %plate
   semicircularbeam layer1 [plate_bottom plate_right]  [radius=radius_plate w=w_plate h=h_plate alpha=arcplate oz=0*oz_plate]
   semicircularbeam layer1 [plate_right  plate_top]    [radius=radius_plate w=w_plate h=h_plate alpha=arcplate oz=1*oz_plate]
   semicircularbeam layer1 [plate_top    plate_left]   [radius=radius_plate w=w_plate h=h_plate alpha=arcplate oz=2*oz_plate]
   semicircularbeam layer1 [plate_left   plate_bottom] [radius=radius_plate w=w_plate h=h_plate alpha=arcplate oz=3*oz_plate]  
   %rim
   semicircularbeam layer1 [plate_bottom plate_right]  [radius=radius_rim w=w_rim h=h_rim alpha=arcplate oz=0*oz_plate oy1=pi/2 oz1=(-pi/2+theta_plate) L1=L1_plate oy2=theta_plate oz2=pi L2=L2_plate] 
   semicircularbeam layer1 [plate_right  plate_top]    [radius=radius_rim w=w_rim h=h_rim alpha=arcplate oz=1*oz_plate oy1=pi/2 oz1=(-pi/2+theta_plate) L1=L1_plate oy2=theta_plate oz2=pi L2=L2_plate]
   semicircularbeam layer1 [plate_top    plate_left]   [radius=radius_rim w=w_rim h=h_rim alpha=arcplate oz=2*oz_plate oy1=pi/2 oz1=(-pi/2+theta_plate) L1=L1_plate oy2=theta_plate oz2=pi L2=L2_plate]
   semicircularbeam layer1 [plate_left   plate_bottom] [radius=radius_rim w=w_rim h=h_rim alpha=arcplate oz=3*oz_plate oy1=pi/2 oz1=(-pi/2+theta_plate) L1=L1_plate oy2=theta_plate oz2=pi L2=L2_plate]
]

plate layer1 [plate_left plate_right plate_top][]

% H I N G E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hinge parameters
w_hinge2plate=20u
w_hinge=9u
h_hinge=highSCS
h_hinge2plate=highSCS
l_hingearm=340u-w_hinge2plate/2-w_hinge/2
l_hingecenter=30u/2
l_hinge2plate=60u
theta_hinge2plate=atan((h_hinge2plate/2-h_plate/2)/(w_plate/2))
L2_hinge2plate=sqrt((w_plate/2)^2+(h_hinge2plate/2-h_plate/2)^2)
%tether parameters
l_tether=735u
h_tether=lowSCS
w_tether=12u

%hinge, left 1 and right 2 
hinge1 layer1 [plate_left hinge1_center hinge1_left][] 
hinge2 layer1 [plate_right hinge_center2 hinge_right2][]
%tether for the hinge, left and right
beam3dlink layer1 [hinge_center2 bridgeright][l=l_tether w=w_tether h=h_tether oz=-pi/2 L1=h_hinge/2-h_tether/2 oy1=pi/2 L2=h_hinge/2-h_tether/2 oy2=-pi/2] 
beam3dlink layer1 [hinge1_center bridgeleft][l=l_tether w=w_tether h=h_tether oz=-pi/2 L1=h_hinge/2-h_tether/2 oy1=pi/2 L2=h_hinge/2-h_tether/2 oy2=-pi/2] 

subnet hinge1 [plate_left hinge1_center hinge1_left][]
[
   %beam from hinge to plate
   beam3dlink layer1 [hinge1_right plate_left][l=l_hinge2plate w=w_hinge2plate h=h_hinge2plate L1=w_hinge/2 L2=L2_hinge2plate oy2=-theta_hinge2plate] %L2 is rigid link to plate   
   %hinge
   beam3dlinkcorner layer1 [hinge1_right  hinge1_rightcorner] [l=l_hingearm    w=w_hinge h=h_hinge oz=pi/2 L1=w_hinge2plate/2 L2=w_hinge/2] %Right arm
   beam3dlink layer1       [hinge1_center hinge1_rightcorner] [l=l_hingecenter w=w_hinge h=h_hinge L2=w_hinge/2] %Center right
   beam3dlink layer1       [hinge1_leftcorner  hinge1_center] [l=l_hingecenter w=w_hinge h=h_hinge L1=w_hinge/2] %Center left
   beam3dlinkcorner layer1 [hinge1_leftcorner    hinge1_left] [l=l_hingearm    w=w_hinge h=h_hinge oz=-pi/2 L1=w_hinge/2 L2=w_hinge2plate/2] %Left arm
]

subnet hinge2 [plate_right hinge_center2 hinge_right2][]
[
   %beam from hinge to plate
   beam3dlink layer1 [hinge_left plate_right][l=l_hinge2plate w=w_hinge2plate h=h_hinge2plate L1=w_hinge/2 L2=L2_hinge2plate oy2=-theta_hinge2plate oz=pi] %L2 is rigid link to plate

   %hinge
   beam3dlinkcorner layer1 [hinge_left  hinge_leftcorner] [l=l_hingearm    w=w_hinge h=h_hinge oz=pi/2 L1=w_hinge2plate/2 L2=w_hinge/2] %Left arm
   beam3dlink layer1 [hinge_leftcorner hinge_center2] [l=l_hingecenter w=w_hinge h=h_hinge L1=w_hinge/2] %Center left
   beam3dlink layer1 [hinge_center2 hinge_rightcorner] [l=l_hingecenter w=w_hinge h=h_hinge L2=w_hinge/2] %Center right
   beam3dlinkcorner layer1 [hinge_rightcorner    hinge_right2] [l=l_hingearm    w=w_hinge h=h_hinge oz=-pi/2 L1=w_hinge/2 L2=w_hinge2plate/2] %Right arm
]

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

subnet perforatedbeam1 [hinge1_left][]
[
   %beam from perforated to hinge
   beam3dlink layer1 [perf_right hinge1_left][l=l_perf2hinge w=w_perf2hinge h=h_perf2hinge L2=w_hinge/2 L1=(h_perf2hinge/2-h_perf/2) oy1=pi/2]

   %short and long perforated beam
   perforated layer1 [perf_middle perf_right][numberholes=perfholes-20 w=w_perf h=h_perf l=l_perf1 whorizontal=w_holehori wvertical=w_holevert] 
   perforated layer1 [perf_left  perf_middle][numberholes=perfholes-4  w=w_perf h=h_perf l=l_perf2 whorizontal=w_holehori wvertical=w_holevert] 
 
   %overhang of perforated
   beam3dlink layer1 [perf_middle perf_right][l=l_perfoverhang h=h_perfoverhang w=w_perfoverhang L1=L_perfoverhang oz1=-theta_perfoverhang oy1=pi/2 L2=L_perfoverhang oz2=theta_perfoverhang oy2=-pi/2] 
   beam3dlink layer1 [perf_middle perf_right][l=l_perfoverhang h=h_perfoverhang w=w_perfoverhang L1=L_perfoverhang oz1=theta_perfoverhang oy1=pi/2 L2=L_perfoverhang oz2=-theta_perfoverhang oy2=-pi/2] 
 
   %anchor for perforated beam
   anchor layer1 [perf_left][l=50u h=upperSCS w=100u oz=pi]  
]

perforatedbeam1 layer1 [hinge1_left][]

% B R I D G E  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_bridge=w_tether/2+l_hingecenter+w_hinge+l_hinge2plate+w_plate
h_bridge=50u
w_bridge=30u
beam3dlink layer1 [bridgeleft  bridgecenter][l=l_bridge w=w_bridge h=h_bridge L1=-w_tether/2] 
beam3dlink layer1 [bridgecenter bridgeright][l=l_bridge w=w_bridge h=h_bridge L2=-w_tether/2] 

% C O M B D R I V E   A R R A Y  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

perfholes=6
w_perf=30u
h_perf=highSCS
w_horizontal=10u
w_vertical=8u
l_perf=100u

w_finger=3.2u
l_finger=30u
gap=3u
numberfingers=70
w_perf=30u
whorizontal=(w_perf-12u)/2
h_perf=highSCS
w_vertical=6u

V=0u
Lcomb=w_bridge/2

perforated layer1 [bridgecenter comb1][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 

subnet combarray [comb1 comb4][]
[
   perforated layer1 [bridgecenter comb1][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 
   perforatedcomb4 layer1 [leftcomb1  comb1][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical L1=0 L2=Lcomb] 
   perforatedcomb4 layer1 [comb1 rightcomb1][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical L1=Lcomb L2=0] 

   l_perf=140u
   perforated layer1     [comb1      comb2][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 
   perforatedcomb4 layer1 [leftcomb2  comb2][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical L1=0 L2=Lcomb] 
   perforatedcomb4 layer1 [comb2 rightcomb2][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical L1=Lcomb L2=0] 

   perforated layer1     [comb2      comb3][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 
   perforatedcomb4 layer1 [leftcomb3  comb3][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical L1=0 L2=Lcomb] 
   perforatedcomb4 layer1 [comb3 rightcomb3][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical L1=Lcomb L2=0] 
   
   perforated layer1     [comb3      comb4][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 
   perforatedcomb4 layer1 [leftcomb4  comb4][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical L1=0 L2=Lcomb] 
   perforatedcomb4 layer1 [comb4 rightcomb4][numberfingers=numberfingers w=w_perf h=h_perf whorizontal=whorizontal wfinger=w_finger lfinger=l_finger gap=gap wvertical=w_vertical L1=Lcomb L2=0] 
]

%combarray layer1 [comb1 comb4][]

l_perf=140u
%perforated layer1 [comb4 comb5][numberholes=perfholes  w=w_perf h=h_perf l=l_perf whorizontal=w_horizontal wvertical=w_vertical oz=-pi/2] 



% P R E D I S P L A C E D   B E A M S  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_predisp=1100u
w_predisp=6.5u
h_predisp=highSCS
y2_disp=-25u
displink=39u

predisplacedbeam4 layer1 [bridgecenter disp1][l=l_predisp w=w_predisp h=h_predisp qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=y2_disp qz2=0 qox2=0 qoy2=0 qoz2=0 L1=displink oz1=-pi/2+pi/8] 
anchor layer1 [disp1][l=50u h=highSCS w=100u ]  
%predisplacedbeam4 layer1 [comb5        disp2][l=l_predisp w=w_predisp h=h_predisp qox1=0 qoy1=0 qoz1=0 qx2=0 qy2=y2_disp qz2=0 qox2=0 qoy2=0 qoz2=0 L1=displink oz1=pi/2-pi/8] 
%anchor layer1 [disp2][l=50u h=highSCS w=100u ]  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second H I N G E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

theta_hinge2plate=atan((h_hinge2plate/2-h_plate/2)/(w_plate/2))
L2_hinge2plate=sqrt((w_plate/2)^2+(h_hinge2plate/2-h_plate/2)^2)

subnet hinge2 [plate_right hinge_center2 hinge_right2][]
[
   %beam from hinge to plate
   beam3dlink layer1 [hinge_left plate_right][l=l_hinge2plate w=w_hinge2plate h=h_hinge2plate L1=w_hinge/2 L2=L2_hinge2plate oy2=-theta_hinge2plate oz=pi] %L2 is rigid link to plate

   %hinge
   beam3dlinkcorner layer1 [hinge_left  hinge_leftcorner] [l=l_hingearm    w=w_hinge h=h_hinge oz=pi/2 L1=w_hinge2plate/2 L2=w_hinge/2] %Left arm
   beam3dlink layer1 [hinge_leftcorner hinge_center2] [l=l_hingecenter w=w_hinge h=h_hinge L1=w_hinge/2] %Center left
   beam3dlink layer1 [hinge_center2 hinge_rightcorner] [l=l_hingecenter w=w_hinge h=h_hinge L2=w_hinge/2] %Center right
   beam3dlinkcorner layer1 [hinge_rightcorner    hinge_right2] [l=l_hingearm    w=w_hinge h=h_hinge oz=-pi/2 L1=w_hinge/2 L2=w_hinge2plate/2] %Right arm
]

hinge2 layer1 [plate_right hinge_center2 hinge_right2][]

%tether for the hinge
beam3dlink layer1 [hinge_center2 bridgeright][l=l_tether w=w_tether h=h_tether oz=-pi/2 L1=h_hinge/2-h_tether/2 oy1=pi/2 L2=h_hinge/2-h_tether/2 oy2=-pi/2] 

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


subnet perforatedbeam2 [hinge_right2][]
[
   %beam from hinge to perforated 
   beam3dlink layer1 [hinge_right2 perf_left][l=l_perf2hinge w=w_perf2hinge h=h_perf2hinge L1=w_hinge/2 L2=(h_perf2hinge/2-h_perf/2) oy2=-pi/2]

   %short and long perforated beam
   perforated layer1 [perf_left perf_middle][numberholes=perfholes-20 w=w_perf h=h_perf l=l_perf1 whorizontal=w_holehori wvertical=w_holevert] 
   perforated layer1 [perf_middle perf_right][numberholes=perfholes-4  w=w_perf h=h_perf l=l_perf2 whorizontal=w_holehori wvertical=w_holevert]    
   
   %overhang 
   beam3dlink layer1 [perf_left perf_middle][l=l_perfoverhang h=h_perfoverhang w=w_perfoverhang L1=L_perfoverhang oz1=-theta_perfoverhang oy1=pi/2 L2=L_perfoverhang oz2=theta_perfoverhang oy2=-pi/2] 
   beam3dlink layer1 [perf_left perf_middle][l=l_perfoverhang h=h_perfoverhang w=w_perfoverhang L1=L_perfoverhang oz1=theta_perfoverhang oy1=pi/2 L2=L_perfoverhang oz2=-theta_perfoverhang oy2=-pi/2] 

   %anchor for perforated beam
   anchor layer1 [perf_right][l=50u h=upperSCS w=100u]  
]

perforatedbeam2 layer1 [hinge_right2][]

