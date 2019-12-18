
%nanosonix netlist
%
%
%
%p.dT=0.01; p.gap=2e-6; p.Nf=10; p.L=200e-6; p.h=2e-6; p.w=2e-6; p.angle=10; p.V=0; net=cho_load('nanosonix_design_BC_02.m',p);figure(1);q=cho_dc(net); cho_display(net,q); [Displacement, Change_in_Capacitance] = nanosonix_design_BC_cap_disp(net,q,p.h,p.gap,p.Nf)

uses mumps.net
%Ambient temperature is 300 Kelvin.

uses subnet_comb_drive_3.m

param V = 10
param w = 2u
param angle = 10 
angle = angle * pi/180
param L = 200u
param dT = 0
T = 300+dT
param gap = 2u
param h = 2u
param nf = 10
wp1 = 2u
wp0 = 0.5u
wsupport = 10u
OL = 5u
Lbackbone = 20u
lf = 20u

%Comb drives
subnet_comb_drive_3 p1 [dn B][oz=0 V=V Nf=nf gapf=gap Lf=lf Wf=w Hf=h overlap=OL W_support=wsupport  ] 
subnet_comb_drive_3 p1 [up C][oz=pi V=V Nf=nf gapf=gap Lf=lf Wf=w Hf=h overlap=OL W_support=wsupport  ] 
%beam3d p1 [dn B] [l=100u w=2u oz=-pi/2]
%beam3d p1 [up C] [l=100u w=2u oz=pi/2]
%Bond pads
bondsize = 100u
bondingpad  p1 [B][l=bondsize w=bondsize h=1.9u oz=-pi/2]
bondingpad  p1 [C][l=bondsize w=bondsize h=1.9u oz=pi/2]

%Backbone
beam3d  p1 [middle up][l=Lbackbone w=10u h=h+1n  oz=pi/2]
beam3d  p1 [middle dn][l=Lbackbone w=10u h=h+1n  oz=-pi/2]

%Flexures 
anchorsize = 20u
beam3d  p1 [left_anchor middle][l=L w=w h=h  oz=angle T=T]
bondingpad  p1 [right_anchor][l=anchorsize w=anchorsize h=4u oz=0]
beam3d  p1 [middle right_anchor ][l=L w=w h=h  oz=-angle T=T]
bondingpad  p1 [left_anchor][l=anchorsize w=anchorsize h=4u oz=pi]

%Oxide
beam3d  p1 [B sb] [h=10u l=2u+wp1/2+wp0/2 w=10u oy=pi/2]

%P0 ground layer
Len2 = wsupport/2 + (lf-OL*2)/2
Len = (wsupport/2 + lf-OL*2 + lf + wsupport + Lbackbone)*2 - Len2*2
beam3dlink  p0 [sb sm] [h=0.5u l=Len w=L*2+anchorsize*2 oz=pi/2 L1=Len2 L2=Len2]

%To tracer
Len3 = (wsupport/2 + lf-OL*2 + lf + wsupport + Lbackbone)
commonground  p0 [right_anchor] [h=0.5u l=Len3-L*sin(angle)+bondsize+10u w=2u oz=-pi/2]
%commonground  p0 [tracer][l=anchorsize w=anchorsize h=4u oz=-pi/2]

