% net=cho_load('circuit1.m');figure(1);cho_display(net);
% net=cho_load('circuit1.m');q=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','y'))

uses process_mumps_cmos1.m

blen = 100u
anchor  p1 [c7]  [l=20u w=20u h=4u]
beam3d p1 [c7 c] [l=blen w=4u]
beam3d p1 [c d]  [l=blen/8 w=4u oz=-pi/2]
beam3d p1 [d d1]  [l=blen/2 w=8u oz=0]

anchor  p1 [c8]  [l=20u w=20u h=4u]
beam3d p1 [c8 c9] [l=blen/2 w=4u]
beam3d p1 [c9 c10] [l=blen/2 w=8u]


len = 30u
resistor  cmos [a b] [len=len R=55 oz=0 L1=len/2 L2=len/2 oz1=-pi/2 oz2=pi/2]

resistor cmos [b c1] [len=len R=1M oz=0 L1=len/2 oz1=pi/2 L2=0u oz2=pi/2]
voltage   cmos [gnd a] [len=len V=15 oz=0 L1=len/2 oz1=pi/2 L2=0u oz2=0]
voltagesinusoidal cmos [aa a2] [len=len V=15 oz=0 L1=len/2 oz1=pi/2 L2=len/2 oz2=pi/2]
resistor cmos [a2 a3] [len=len R=55 oz=0 L1=len/2 L2=len/2 oz1=0 oz2=0]

capacitor cmos [c1 c6] [len=len C=55p oz=0 L1=0u L2=len/2 oz1=0 oz2=pi/2]
inductor cmos [a3 a4] [len=len L=33 oz=-pi/2 L1=0*len oz1=0 L2=0*len oz2=-pi/2]
ground cmos [gnd eg] [len=len V=0 oz=-pi/2]
short cmos [aa a5] [len=2*len V=0 oz=0 L2=len/2 oz2=-pi/2]
short cmos [b c7] [len=len/2 V=0 oz=0 L2=len/2 oz2=-pi/4]
eswitch cmos [a aa] [len=len oz=pi/2]

short cmos [gnd c8] [len=len V=0 oz=-pi/4 L1=len/2 oz1=pi/4 L2=99.5u oz2=pi/4]

