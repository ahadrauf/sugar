% net=cho_load('circuit3.m');figure(1);cho_display(net);
% net=cho_load('circuit3.m');q=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','y'))

uses process_mumps_cmos1.m
blen = 100u
len = 30u

anchor     p1 [c70]   [l=40u w=40u h=4u oz=pi]
beam3dlink p1 [c70 c] [l=200u w=4u L1=18u oz1=-pi/2 L2=2u oz2=-pi/2] %hot wire
beam3dlink p1 [c d]   [l=4u w=4u oz=-pi/2 L1=2u oz1=-pi/2 L2=2u oz2=pi/2] %link
beam3dlink p1 [d1 d]  [l=160u w=20u L2=10u oz2=pi/2] %cold wire
beam3dlink p1 [d20 d1][l=40u w=4u L2=8u oz2=-pi/2 L1=18u oz1=pi/2] %attachment
anchor     p1 [d20]   [l=40u w=40u h=4u oz=pi]
beam3d     p1 [c7 c70][l=40u w=2u h=1u] 
beam3d     p1 [d2 d20][l=40u w=2u h=1u] 
beam3dlink p1 [d3 d20][l=200u w=0.1u h=0.1u] 

resistor    cmos [a2 c7] [len=len R=55 oz=-pi/2 L2=len oz2=pi/2]
capacitor   cmos [a3 a2] [len=len C=55p oz=0 L2=len/2 oz2=-pi/2]
inductor    cmos [a4 a3] [len=len L=33 L2=len/4 L1=len/4]
voltagesinusoidal cmos [a5 a4] [len=len V=15 ]
resistor    cmos [a4 a6] [len=len V=0 oz=-pi/2 L1=len/4 L2=len/4]
eswitch     cmos [a6 a7] [len=len oz=0]
short       cmos [a7 a8] [len=3*len/4 V=0 oz=-pi/2 L1=len/2 oz1=pi/2 L2=len/2 oz2=-pi/2]
voltage     cmos [b1 a8] [len=len V=15 oz=0 L1=len]

resistor    cmos [a9 d2] [len=len R=55 oz=0 L2=len/2 L1=len/2]
capacitor   cmos [a9 a10][len=len C=55p oz=-pi/2]
ground      cmos [a10 gnd1] [len=len V=0 oz=-pi/2]
short       cmos [a11 a9] [len=len*1.8 V=0 oz=0 ]

short       cmos [a5 b1] [len=len*2.25 V=0 oz=-pi/2]
ground      cmos [a11 gnd2] [len=len V=0 oz=-pi/2 L1=len]
short       cmos [b1 a11] [len=len V=0 oz=-pi/4]

