%net=cho_load('robotleg2.m');figure(2);cho_display(net);shownodes(net);view(2)
%net=cho_load('robotleg2.m');figure(2);q=cho_dc(net);cho_display(net,q);shownodes(net);view(2)

%net=cho_load('robotleg1.m');figure(1);cho_display(net);shownodes(net);view([-1 1 -1])
%net=cho_load('robotleg1.m');figure(1);cho_display(net);
%beam and slider
% net=cho_load('slider1.m');[q,K,Tv,F,T,v]=cho_dc_slider(net);figure(1);cho_display(net,q);
% net=cho_load('slider1.m');figure(1);cho_display(net);
%for f=0:-5e-6:-40e-6, p.f=f;net=cho_load('slider1.m',p);[q]=cho_dc_slider(net);figure(1);cho_display(net,q); end

uses mumps.net
%param f=-30u
%a=pi/3
%anchor p1 [C]   [l=10u w=10u oz=-pi/2] %lower anchor
%beam3d p1 [A C] [l=200u w=2u h=2u oz=-a] %lower beam
%slider2 p1 [A B] [l=200u w=2u h=2u oz=a slideroz=a] %slider on node2
%f3d * [A][F=f] %horizontal force

L1=150u
L2=50u
L3=100u
L5=100u
L6=300u
L7=L2/2
w1=50u
w3=100u
L2=w3/2
w4=L1*0.9
w5=L3+w1
h=10u
wh=1u
hh=1u
offset=6u
tetherprojection=400u
theta=pi/3
alpha=-pi/3
L9=500u
L8=20u
h9=L8*2* 0.9
w9=400u


subnet leg [T1 T2 anch][]
[

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
beam3d p1 [B3 o] [l=10u w=10u h=h ox=pi/2 oz=alpha] 

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
beam3d p1 [b a] [l=L9 w=w9 h=h9 ox=pi/2] 

]

leg p1 [T1 T2 A1][ oz=0 ox=-pi/2]

f3d * [T1][F=-4u] 
f3d * [T2][F=2u] 

