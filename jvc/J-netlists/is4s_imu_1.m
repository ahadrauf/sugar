% net=cho_load('is4s_imu_1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','z'))
% net=cho_load('is4s_imu_1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('is4s_imu_1.m');figure(1);cho_display(net);

uses mumps.net
mL = 500u
H = 8u
W = 2.11u
W2 = 20u
fL = 200u
aL = 50u 

%anchor p1 [m]   [l=10u w=10u oz=pi]

%proof mass
beam3dlink p1 [m m(1)] [L1=mL/4 oz1=pi/2  l=mL/2 w=mL/2 h=H L2=mL/4 oz2=pi/2] %upper right
beam3dlink p1 [m(4) m] [L2=mL/4 oz2=-pi/2 l=mL/2 w=mL/2 h=H L1=mL/4 oz1=-pi/2] %upper right
beam3dlink p1 [m(7) m] [L2=mL/4 oz2=pi/2  l=mL/2 w=mL/2 h=H L1=mL/4 oz1=pi/2] %lower right
beam3dlink p1 [m m(10)][L1=mL/4 oz1=-pi/2 l=mL/2 w=mL/2 h=H L2=mL/4 oz2=-pi/2] %lower left

%flexures from mass
for i = 1 : 4
[
	beam3d p1 [m(3*i-2) m(3*i-2+2)][l=fL w=W h=H oz=i*pi/2] %1st vertical
    if i == 1
	[
        beam3d p1 [m(3*i-2) m(14)][l=fL w=W h=H oz=(i-1)*pi/2] %1st horizontal
    ]
    else
    [
    	beam3d p1 [m(3*i-2) m(3*i-2+1)][l=fL w=W h=H oz=(i-1)*pi/2] %1st horizontal
    ]
]

%comb drive backbone with flexures
for i = 1 : 4
[
    beam3dlink p1 [m(3*i) c(3*i)][l=mL/2 w=W2 h=H oz=pi+(i-1)*pi/2 oz2=-pi/2 L2=W2/2] %backbone, 1st horizontal, r2l    
    beam3dlink p1 [c(3*i) m(3*i+2)][oz1=pi/2 L1=W2/2 l=mL/2 w=W2 h=H oz=pi+(i-1)*pi/2] %backbone
    beam3d p1 [m(3*i) m1(3*i)][l=fL-W2 w=W h=H oz=0+(i-1)*pi/2] %right flexure
    beam3d p1 [m1(3*i) m2(3*i)][l=W2 w=W h=H oz=pi/2+(i-1)*pi/2] %crab
     anchor p1 [m2(3*i)][l=aL w=aL h=H oz=pi/2+(i-1)*pi/2] %anchor
    beam3d p1 [m(3*i+2) m1(3*i+2)][l=fL-W2 w=W h=H oz=pi+(i-1)*pi/2] %left flexure
    beam3d p1 [m1(3*i+2) m2(3*i+2)][l=W2 w=W h=H oz=pi/2+(i-1)*pi/2] %crab
     anchor p1 [m2(3*i+2)][l=aL w=aL h=H oz=pi/2+(i-1)*pi/2] %anchor
] 
%more comb drive backbones
for i = 1 : 4*4
[
    beam3d p1 [c(3*i) c(3*i+12)][l=W2*3 w=W2 h=H oz=pi/2+(i-1)*pi/2] %coupling  
    beam3dlink p1 [cc(3*i) c(3*i+12)][l=mL/2 w=W2 h=H oz=pi+(i-1)*pi/2 oz2=-pi/2 L2=W2/2] %backbone, 1st horizontal, r2l    
    beam3dlink p1 [c(3*i+12) cc(3*i+2)][oz1=pi/2 L1=W2/2 l=mL/2 w=W2 h=H oz=pi+(i-1)*pi/2] %backbone
]
%Fingers
for j = 1 : 4*5
[    
    for i = 1 : floor((fL/2)/(2*W))
    [
        beam3dlink p1 [c(3*j) c1(j,i)][L1=mL/2*(i/floor((fL/2)/(2*W))) oz1=-pi/2 l=W2 w=W h=H oz=pi/2+(j-1)*pi/2] %finger
        beam3dlink p1 [c(3*j) c2(j,i)][L1=-mL/2*(i/floor((fL/2)/(2*W))) oz1=-pi/2 l=W2 w=W h=H oz=pi/2+(j-1)*pi/2] %finger
    ]
]



f3d * [c(3)][F=65e-6 oz=pi/2] 
f3d * [c(9)][F=65e-6 oz=pi/2] 

%anchor p1 [mf(4)][l=aL w=aL h=H oz=0] %anchor
