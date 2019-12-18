%net=cho_load('robotleg4.m');figure(1);cho_display(net);

uses mumps.net

thighfactor=1.5
shinfactor=1.5
L1=150u*thighfactor
L2=50u*thighfactor
L3=100u
L5=100u*shinfactor
L6=300u*shinfactor
L7=L2/2
w1=50u
w3=L2*2*0.8
L2=w3/2
w4=L1*0.9
w5=L3+w1
h=20u
wh=1u
hh=1u
offset=6u
tetherprojection=400u
theta=pi/3
alpha=-pi/3
L9=500u
L8=20u
h9=L8*2* 0.9
w9=400u*1.3

subnet leg [T1 T2 R][fT1=* fT2=*]
[
   
f3d * [T1][F=fT1] 
f3d * [T2][F=fT2] 

%body hinge
anchor p1 [anch]   [l=L9 w=w9 h=5u oz=pi ox=pi/2] 
beam3d p1 [ALl anch] [l=L3/2 w=wh h=hh oy=-pi/2] 
beam3d p1 [anch ALr] [l=L3/2 w=wh h=hh oy=-pi/2] 
%body legs attached to hinge
beam3d p1 [ALl AMl] [l=L1 w=w1 h=h oz=theta ox=pi/2] 
beam3d p1 [ALr AMr] [l=L1 w=w1 h=h oz=theta ox=pi/2] 
%middle leg hinge
beam3d p1 [AMl h1] [l=L3/2-offset w=wh h=hh oy=-pi/2] 
beam3d p1 [h1 AMr] [l=L3/2+offset w=wh h=hh oy=-pi/2] 
%leg going to knee
beam3d p1 [AMl AHl] [l=L2 w=w1 h=h oz=theta ox=pi/2 ] 
beam3d p1 [AMr AHr] [l=L2 w=w1 h=h oz=theta ox=pi/2 ] 
beam3d p1 [AHl AHr] [l=L3 w=w3*0.8 h=h*0.8 oy=-pi/2 ox=pi/2+theta] %mid section
beam3d p1 [AHl ahl] [l=L2 w=w1 h=h oz=theta ox=pi/2 ] 
beam3d p1 [AHr ahr] [l=L2 w=w1 h=h oz=theta ox=pi/2 ] 
%knee hinge
beam3dlink p1 [ahl h2] [l=L3/2 w=wh h=hh oy=-pi/2 ] 
beam3dlink p1 [h2 Ahr] [l=L3/2 w=wh h=hh oy=-pi/2 ] 

%shin connected to knee
beam3d p1 [h2 B1] [l=L5 w=w5 h=h ox=pi/2 oz=alpha] 
%hinge supports
Llink=w5/2-w1/2
beam3dlink p1 [B1 B4] [l=L7 w=w1 h=h ox=pi/2 L1=Llink oz1=pi/2 oz=alpha] 
beam3dlink p1 [B1 B5] [l=L7 w=w1 h=h ox=pi/2 L1=Llink oz1=-pi/2 oz=alpha] 
beam3dlink p1 [B4 B2] [l=L7 w=w1 h=h ox=pi/2 L2=Llink oz2=-pi/2 oz=alpha] 
beam3dlink p1 [B5 B2] [l=L7 w=w1 h=h ox=pi/2 L2=Llink oz2=pi/2 oz=alpha] 
%shin hinge
beam3d p1 [B5 h3] [l=L3/2+offset w=wh h=hh oy=-pi/2] 
beam3d p1 [h3 B4] [l=L3/2-offset w=wh h=hh oy=-pi/2] 
%foot
beam3d p1 [B2 B3] [l=L6 w=w5 h=h ox=pi/2 oz=alpha] 
beam3d p1 [B3 o] [l=40u w=40u h=h ox=pi/2 oz=alpha] 

%tethers
Tv1=L1*sin(theta)
Th1=tetherprojection
Tv2=(L1+2*L2)*sin(theta)+(L5+L7)*sin(alpha)
Th2=tetherprojection
tw=3u
th=3u
slider3 p1 [h1 T1] [l=sqrt(Tv1^2+Th1^2) w=tw h=th oz=pi+atan(Tv1/Th1) ] %slider on node 2
slider3 p1 [h3 T2] [l=sqrt(Tv2^2+Th2^2) w=tw h=th oz=pi+atan(Tv2/Th2) ] 
%slider attachments
slider3 p1 [T1 TT1] [l=150u w=40u h=5u oz=pi]
slider3 p1 [T2 TT2] [l=150u w=40u h=5u oz=pi]

%body
beam3d p1 [anch a] [l=L8 w=1u h=1u oz=-pi/2] 
beam3d p1 [R a] [l=L9 w=w9 h=h9 ox=pi/2] 

]

%mandible right
m1=20u
hm=30u
lm=100u
lmh=80u
subnet hand [m1][]
[
   beam3d p1 [m1 m2][l=lm w=m1 h=hm oz=(pi/2+pi/6)]
   beam3d p1 [m4 m2][l=lmh w=m1 h=hm ]
   beam3d p1 [m5 m3][l=lmh w=m1 h=hm ]
   beam3d p1 [m2 m3][l=lmh w=m1 h=hm oz=pi/2]
]
subnet mandible [R1][side=* oyh=*]
[
   %arm
   beam3dlink p1 [R1 m1][l=lm w=m1 h=hm oz=0 L1=sqrt((w9/2)^2+(L9/3)^2) oz1=side*atan(L9/w9) oz=pi/2]
   hand p1 [m1][oy=oyh]
]


subnet bot [a] []
[
ox=0
oz=0
oy=0
% fT1 => thigh 
% fT2 => shin 
% negative force => retract
% positive force => extend
leg p1 [T11b T12b R1][ oz=oz ox=ox oy=oy			fT1=4u fT2=0.5u] 
%leg p1 [T11a T12a R1][ oz=oz ox=ox oy=oy+pi		fT1=0u fT2=0u]
%beam3d p1 [R1 R2] [l=w9 w=1u h=1u oy=pi/2] 
%leg p1 [T21b T22b R2][ oz=oz ox=ox oy=oy			fT1=0u fT2=0u]
%leg p1 [T21a T22a R2][ oz=oz ox=ox oy=oy+pi		fT1=2u fT2=2u]
%beam3d p1 [R2 R3] [l=w9 w=1u h=1u oy=pi/2] 
%leg p1 [T31b T32b R3][ oz=oz ox=ox oy=oy			fT1=0u fT2=0u]
%leg p1 [T31a T32a R3][ oz=oz ox=ox oy=oy+pi		fT1=0u fT2=0u]
%mandible
%mandible p1 [R1][ox=pi/2 side=-1 ox=pi/2 oyh=0]
%mandible p1 [R1][ox=pi/2 side=1 ox=pi/2 oyh=pi]
]

bot p1 [a][ox=pi]
