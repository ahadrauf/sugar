% net=cho_load('nanosonix_DesignA_3.net',p); q=cho_dc(net); figure(1); clf; cho_display(net,q); 
%one = q(lookup_coord(net, 'c', 'y'))

uses mumps.net

%Parameters to modify
L = 100u
w = 10u

%Other parameters
L1 = 4u
hp1 = 2u
ho2 = 2u
hp0 = 0.5u
Lanchor = 15u
wanchor = Lanchor
Dw=10u
welec=w+2*Dw
tipgap = 2u

%anchor 
bondingpad  p1 [a] [l=Lanchor w=wanchor oz=-pi]
beam3dlink_cloak p1 [a c] [L1=Lanchor/2 l=(hp1/2 + ho2 + hp0) w=wanchor h=wanchor oy1=pi/2 oy=pi/2]

%Cantilever
beam3d p1 [a b] [l=L w=w ]
f3d * [b] [ M=1p oz=pi/2 ]

%Electrodes
L5=(2*L/4-L1)/2
L3=(L/4-L1)/2
L8=L1
beam3dlink p0 [c d] [L1=L1+Lanchor/2 l=L3*2 w=welec ]
beam3dlink p0 [d e] [L1=L1 l=2*L5 w=welec ]
beam3dlink p0 [e f] [L1=L1 l=L3*2 w=welec L2=L1]

%Electrode bondpads
bondsize = 100u
L7=Dw
L4=bondsize/2
L6=(welec/2 + Dw)
L9=L4+L7+L4-L3-L8-L5

beam3dlink p0 [d t1] [L1=L3 oz1=pi/2 oz=pi/2 l=L6+Dw w=2u ]
bondingpad  p1 [t1][l=bondsize w=bondsize h=1.9u oz=pi/2]

beam3dlink p0 [e e2] [L1=L5 oz1=pi/2 oz=pi/2 l=L6 w=2u ]
   beam3dlink p0 [e2 ee2] [oz1=pi/2 oz=pi/2 l=1u w=2u ]
beam3dlink p0 [e2 e3] [L1=1u l=L9-1u w=2u ]
beam3dlink p0 [e3 t2] [L1=-1u l=Dw+1u w=2u oz=pi/2]
bondingpad  p1 [t2][l=bondsize w=bondsize h=1.9u oz=pi/2]

L10=L3+L8+L5*2+L8+L3
beam3dlink p0 [f f2] [L1=L3+L1 oz1=-pi/2 oz=-pi/2 l=L6 w=2u ]
beam3dlink p0 [f2 f3] [L1=-1u l=L10 w=2u oz=pi]
beam3dlink p0 [f3 t3] [L1=-1u l=Dw+1u w=2u oz=-pi/2]
bondingpad  p1 [t3][l=bondsize w=bondsize h=1.9u oz=-pi/2]

L12=L+tipgap+L/10+Lanchor/2
L13=Lanchor/2
L14=Lanchor/2+L/10+tipgap+L3*2+L8+L5-L9
beam3dlink p0 [a a4] [L1=L12 oz1=pi/2 oz=-pi/2 l=L6+Dw+2u w=2u L2=-1u]
beam3dlink p0 [t4 a4] [l=L14 w=2u L1=1u oz1=-pi/2]
bondingpad  p1 [t4][l=bondsize w=bondsize h=1.9u oz=-pi/2]

%Tip electrode
beam3dlink p1 [a g] [L1=L+tipgap  l=L/10 w=w ]
bondingpad  p1 [g] [l=Lanchor w=Lanchor ]
beam3dlink_cloak p1 [g h] [L1=Lanchor/2 l=(hp1/2 + ho2 + hp0) w=wanchor h=wanchor oy1=-pi/2 oy=pi/2]

%To tracer
beam3dlink p0 [T c] [L1=-1u l=bondsize/2-Lanchor/2-Dw w=2u L2=Lanchor/2 ]
commonground  p0 [T gnd] [l=L6+Dw+bondsize+Dw*2 w=2u oz=-pi/2]

%SugarCube parameters
sugarcube * [] [ V = 'Voltage, V, 10, 0, 20' 
	w = 'Cantilever width, m, 10e-6, 2e-6, 50e-6' 
    L = 'Cantilever length, m, 100e-6, 30e-6, 500e-6' 
	nodes = 'b' ]

                        
