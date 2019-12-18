%net=cho_load('robot2.m');q=cho_dc(net);figure(1);cho_display(net,q);  
%net=cho_load('robot2.m');figure(1);cho_display(net);  
 
uses mumps2.net

L1=100u %lenght of plates
W1=200u %width of plates
W2=100u %width of bottom springs
W3=2u %width of top springs
L4=W1/2 %length of springs
g1=100u %leg1 length
g2=400u %leg2 length
gy=60 %leg angles
hp=20u %plate thickness
hl=20u %leg thickness
wl=40u %leg width
pi=3.141592653589793
theta=gy
%T=sqrt(L1^2 + g1^2 + 2*L1*g1*cos(theta*pi/180)) %length of tether
T=1.732050807568877e-004 %
%alpha=asin(g1/T * sin(pi - theta*(pi/180))) 
alpha=30 %tether angle

%right 
beam3d p1 [a  ar] [l=L1 w=W1 h=hp oz=0]
freehinge p1 [ar arr][l=L1 w=W1 h=hp oz=0]
beam3d p1 [b  br] [l=L1 w=W1 h=hp oz=0]
beam3d p1 [br brr][l=L1 w=W1 h=hp oz=0]
beam3d p1 [c  cr] [l=L1 w=W1 h=hp oz=0]
freehinge p1 [cr crr][l=L1 w=W1 h=hp oz=0]
%left
beam3d p1 [a  al] [l=L1 w=W1 h=hp oz=180]
beam3d p1 [al all][l=L1 w=W1 h=hp oz=180]
beam3d p1 [b  bl] [l=L1 w=W1 h=hp oz=180]
freehinge p1 [bl bll][l=L1 w=W1 h=hp oz=180]
beam3d p1 [c  cl] [l=L1 w=W1 h=hp oz=180]
beam3d p1 [cl cll][l=L1 w=W1 h=hp oz=180]

%body weight
beam3d p1 [b bd][l=hp/2*4 w=2u oy=-90]
beam3d p1 [bd bd1][l=W1*1.5 w=L1*4 h=hp*4 oz=90]
beam3d p1 [bd bd2][l=W1*1.5 w=L1*4 h=hp*4 oz=-90]

%right springs 1
beam3d p1 [ar r1] [l=L4 w=W2 oz=90]
beam3d p1 [ar r3] [l=L4 w=W2 oz=-90]
beam3d p1 [r2 r1] [l=L4 w=W3 oz=90]
beam3d p1 [r2 r3] [l=L4 w=W3 oz=-90]
%right springs 2
beam3d p1 [br r3] [l=L4 w=W2 oz=90]
beam3d p1 [br r5] [l=L4 w=W2 oz=-90]
beam3d p1 [r4 r3] [l=L4 w=W3 oz=90]
beam3d p1 [r4 r5] [l=L4 w=W3 oz=-90]
%right springs 3
beam3d p1 [cr r5] [l=L4 w=W2 oz=90]
beam3d p1 [cr r7] [l=L4 w=W2 oz=-90]
beam3d p1 [r6 r5] [l=L4 w=W3 oz=90]
beam3d p1 [r6 r7] [l=L4 w=W3 oz=-90]

%left springs 1
beam3d p1 [al l1] [l=L4 w=W2 oz=90]
beam3d p1 [al l3] [l=L4 w=W2 oz=-90]
beam3d p1 [l2 l1] [l=L4 w=W3 oz=90]
beam3d p1 [l2 l3] [l=L4 w=W3 oz=-90]
%left springs 2
beam3d p1 [bl l3] [l=L4 w=W2 oz=90]
beam3d p1 [bl l5] [l=L4 w=W2 oz=-90]
beam3d p1 [l4 l3] [l=L4 w=W3 oz=90]
beam3d p1 [l4 l5] [l=L4 w=W3 oz=-90]
%left springs 3
beam3d p1 [cl l5] [l=L4 w=W2 oz=90]
beam3d p1 [cl l7] [l=L4 w=W2 oz=-90]
beam3d p1 [l6 l5] [l=L4 w=W3 oz=90]
beam3d p1 [l6 l7] [l=L4 w=W3 oz=-90]

%right leg 1
beam3d p1 [arr    rleg11] [l=g1 w=wl h=hl oy=gy]
beam3d p1 [rleg11 rleg12] [l=g1 w=wl h=hl oy=gy]
beam3d p1 [rleg12 rleg13] [l=g2 w=wl h=hl oy=-gy]
%right leg 2
beam3d p1 [brr    rleg21] [l=g1 w=wl h=hl oy=gy]
beam3d p1 [rleg21 rleg22] [l=g1 w=wl h=hl oy=gy]
beam3d p1 [rleg22 rleg23] [l=g2 w=wl h=hl oy=-gy]
%right leg 3
beam3d p1 [crr    rleg31] [l=g1 w=wl h=hl oy=gy]
beam3d p1 [rleg31 rleg32] [l=g1 w=wl h=hl oy=gy]
beam3d p1 [rleg32 rleg33] [l=g2 w=wl h=hl oy=-gy]
%feet
hinge p1 [rleg23][l=20u w=10u h=10u oy=-90]
beam3d p1 [rleg33 rleg33x][l=20u w=10u h=10u oy=-90]
beam3d p1 [rleg13 rleg13x][l=20u w=10u h=10u oy=-90]

%left leg 1
beam3d p1 [all    lleg11] [l=g1 w=wl h=hl oy=180-gy]
beam3d p1 [lleg11 lleg12] [l=g1 w=wl h=hl oy=180-gy]
beam3d p1 [lleg12 lleg13] [l=g2 w=wl h=hl oy=180+gy]
%left leg 2
beam3d p1 [bll    lleg21] [l=g1 w=wl h=hl oy=180-gy]
beam3d p1 [lleg21 lleg22] [l=g1 w=wl h=hl oy=180-gy]
beam3d p1 [lleg22 lleg23] [l=g2 w=wl h=hl oy=180+gy]
%left leg 3
beam3d p1 [cll    lleg31] [l=g1 w=wl h=hl oy=180-gy]
beam3d p1 [lleg31 lleg32] [l=g1 w=wl h=hl oy=180-gy]
beam3d p1 [lleg32 lleg33] [l=g2 w=wl h=hl oy=180+gy]
%feet
hinge p1 [lleg13][l=20u w=10u h=10u oy=-90]
hinge p1 [lleg33][l=20u w=10u h=10u oy=-90]
beam3d p1 [lleg23 lleg23x][l=20u w=10u h=10u oy=-90]

%tethers
beam3d p1 [r2 rleg11][l=T w=2u h=2u oy=alpha] 
beam3d p1 [r4 rleg21][l=T w=2u h=2u oy=alpha] 
beam3d p1 [r6 rleg31][l=T w=2u h=2u oy=alpha] 
%tethers
beam3d p1 [l2 lleg11][l=T w=2u h=2u oy=180-alpha] 
beam3d p1 [l4 lleg21][l=T w=2u h=2u oy=180-alpha] 
beam3d p1 [l6 lleg31][l=T w=2u h=2u oy=180-alpha] 

%mandible right
wa=10u
ha=20u
la=20u
beam3d p1 [r1 o1]  [l=la*7 w=wa h=ha oz=90]
beam3d p1 [o1 o2]  [l=la*3 w=wa h=ha oz=90+20]
beam3d p1 [o2 o3]  [l=la w=wa h=ha oz=90]
beam3d p1 [o2 o4]  [l=la w=wa h=ha oz=-90]
beam3d p1 [o3 o33] [l=la w=wa h=ha oz=180]
beam3d p1 [o4 o44] [l=la w=wa h=ha oz=180]
%mandible left
beam3d p1 [l1 o21]  [l=la*7 w=wa h=ha oz=90]
beam3d p1 [o21 o22]  [l=la*3 w=wa h=ha oz=90-20]
beam3d p1 [o22 o23]  [l=la w=wa h=ha oz=90]
beam3d p1 [o22 o24]  [l=la w=wa h=ha oz=-90]
beam3d p1 [o23 o233] [l=la w=wa h=ha oz=0]
beam3d p1 [o24 o244] [l=la w=wa h=ha oz=0]
%torso forces
f=100u
f3d * [r2] [F=-f oz=0]
f3d * [l4] [F=f*2 oz=0]

