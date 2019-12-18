% net=cho_load('circuit2.m');figure(1);cho_display(net);
% net=cho_load('circuit2.m');q=cho_dc(net);figure(1);cho_display(net,q);q(lookup_coord(net,'B','y'))

uses process_mumps_cmos1.m
len=30u
%resistor  cmos [a b] [len=len R=55 oz=pi/4 linewidth=5 nodewidth=15]
%capacitor cmos [b c1] [len=len C=44 oz=0 linewidth=5 nodewidth=15]
%voltage   cmos [gnd a] [len=len V=15 oz=0 linewidth=5 nodewidth=15]
%voltagesinusoidal cmos [aa a2] [len=len V=15 oz=0 L1=len/2 oz1=pi/2 L2=len/2 oz2=pi/2 linewidth=5 nodewidth=15]
%inductor cmos [a3 a4] [len=len L=33 oz=0 L1=0*len oz1=0 L2=0*len oz2=-pi/2 linewidth=5 nodewidth=15]
%ground cmos [gnd eg] [len=len V=0 oz=-pi/4 linewidth=5 nodewidth=15]
%short cmos [aa a5] [len=2*len V=0 oz=0 L2=len/2 oz2=-pi/2]
eswitch cmos [a aa] [len=len oz=0 linewidth=5 nodewidth=15]
