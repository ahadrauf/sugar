%micromirror
%net=cho_load('heatuator1.m');figure(1);cho_display(net);
%net=cho_load('heatuator1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%clear all;tic;net=cho_load('heatuator1.m');figure(1);cho_display(net);toc
%clear all;tic;net=cho_load('heatuator1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);toc

% clear all;tic;net=cho_load('heatuator1.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);toc,theta=q(lookup_coord(net,'plate_top','rx')),cosine=q(lookup_coord(net,'bridgecenter','y'))

%tic;net=cho_load('heatuator1.m');toc,[f,egv,dq]=cho_mode(net);toc,mode=1;cho_modeshape(net, f, egv, dq, 5,mode);toc

uses mumps.net

l_anchor=10u
w_anchor=10u
h_anchor=2u

l_hot=150u
w_hot=2u
h_hot=h_anchor

l_bridge=4u
w_bridge=2u
h_bridge=h_anchor

l_cold=l_hot*0.75
w_cold=20u
h_cold=h_anchor

l_short=l_hot*0.25
w_short=2u
h_short=h_anchor



anchor p1 [A] [l=l_anchor w=w_anchor h=h_anchor oz=pi]
anchor p1 [B] [l=l_anchor w=w_anchor h=h_anchor oz=pi]
beam3dlink p1 [A a] [l=l_hot w=w_hot h=h_hot L1=w_anchor/2-w_hot/2 L2=-w_bridge/2 oz1=-pi/2 ] %hot arm   
beam3dlink p1 [a c] [l=l_bridge w=w_bridge h=h_bridge oz=-pi/2 L1=w_hot/2 L2=w_cold/2] %bridge
beam3dlink p1 [b c] [l=l_cold w=w_cold h=h_cold L2=-w_bridge/2] %cold arm
beam3dlink p1 [B b] [l=l_short w=w_short h=h_short L1=w_anchor/2-w_short/2 oz1=pi/2 L2=w_cold/2-w_short/2 oz2=-pi/2 ] %short arm

%f3d * [c] [F=70u oz=pi/2]
f3d * [a] [F=-5000u ]

