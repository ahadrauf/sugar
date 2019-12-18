
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


anchor p1 [ALl]   [l=20u w=20u h=2u oz=pi ox=pi/2] 

L1=150u
L2=50u
L3=100u
w1=50u
w3=50u
w4=L1*0.9
h=10u
%2 low legs
beam3d p1 [ALl AMl] [l=L1 w=w1 h=h oz=pi/2 ox=pi/2] 
beam3d p1 [ALr AMr] [l=L1 w=w1 h=h oz=pi/2 ox=pi/2] 

%hinge
beam3d p1 [AMl h1] [l=L3/2+w1/2-6u w=1u h=1u oy=-pi/2] 
beam3d p1 [h1 AMr] [l=L3/2+w1/2+6u w=1u h=1u oy=-pi/2] 
%upper leg
beam3dlink p1 [AMl AHl] [l=L2 w=w1 h=h oz=pi/2 ox=pi/2 L2=w3/2] 
beam3dlink p1 [AMr AHr] [l=L2 w=w1 h=h oz=pi/2 ox=pi/2 L2=w3/2] 
beam3dlinkcorner p1 [AHl AHr]  [l=L3 w=w1 h=h oy=-pi/2 L1=w1/2 L2=w1/2] 
%knee
beam3dlink p1 [AHl BL] [l=L3/2+w1/2 w=1u h=1u oy=-pi/2 L1=w3/2 oz1=pi/2] 
beam3dlink p1 [BL AHr] [l=L3/2+w1/2 w=1u h=1u oy=-pi/2 L2=w3/2 oz2=-pi/2] 


L5=200u
L6=200u
w5=L3+w1*2
w7=w1
L7=L2
%leg
beam3d p1 [BL B1] [l=L5 w=w5 h=h ox=pi/2] 
%hinge supports
L=w5/2-w7/2
beam3dlink p1 [B1 B4] [l=L7 w=w1 h=h ox=pi/2 L1=L oz1=pi/2] 
beam3dlink p1 [B1 B5] [l=L7 w=w1 h=h ox=pi/2 L1=L oz1=-pi/2] 
beam3dlink p1 [B4 B2] [l=L7 w=w1 h=h ox=pi/2 L2=L oz2=-pi/2] 
beam3dlink p1 [B5 B2] [l=L7 w=w1 h=h ox=pi/2 L2=L oz2=pi/2] 
%leg
beam3d p1 [B2 B3] [l=L6 w=w5 h=h ox=pi/2] 
beam3d p1 [B3 o] [l=10u w=10u h=h ox=pi/2] 
%hinge
beam3d p1 [B5 h2] [l=L3/2+w1/2+6u w=1u h=1u oy=-pi/2] 
beam3d p1 [h2 B4] [l=L3/2+w1/2-6u w=1u h=1u oy=-pi/2] 


%tethers
Tv1=L1
Th1=400u
Tv2=L1+L2+w3
Th2=400u
slider1 p1 [h1 T1] [l=sqrt(Tv1^2+Th1^2) w=3u h=3u oz=pi+atan(Tv1/Th1) slideroz=0] 
slider1 p1 [h2 T2] [l=sqrt(Tv2^2+Th2^2) w=3u h=3u oz=pi+atan(Tv2/Th2) slideroz=pi/4] 

f3d * [T2][F=0.1u] 
