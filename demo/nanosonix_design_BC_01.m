%p.dT=0.01; p.gap=2e-6; p.Nf=10; p.L=200e-6; p.h=2e-6; p.w=2e-6; p.angle=10; p.V=0; net=cho_load('nanosonix_design_BC_01.m',p);figure(1);q=cho_dc(net); cho_display(net,q); [Displacement, Change_in_Capacitance] = nanosonix_design_BC_cap_disp(net,q,p.h,p.gap,p.Nf)

uses mumps.net
%Ambient temperature is 300 Kelvin.

uses subnet_comb_drive_3.m

param V = 10
param w = 2u
param angle = 0
angle = angle * pi/180
param L = 200u
param dT = 0
T = 300+dT
param gap = 2u
param h = 2u
param Nf = 10

param p1_pad_l=100u
param p1_pad_w=100u

subnet_comb_drive_3 p1 [dn B][ V=V Nf=Nf gapf=gap Lf=20u Wf=w Hf=h overlap=5u W_support=10u ] 
% bondingpad  p1 [B][l=p1_pad_l w=p1_pad_l h=1.9u oz=-pi/2]
% subnet_comb_drive_3 p1 [up C][ oz=pi V=V Nf=Nf gapf=gap Lf=20u Wf=w Hf=h overlap=5u W_support=10u ] 
% bondingpad  p1 [C][l=p1_pad_l w=p1_pad_l h=1.9u oz=pi/2]

%Anchor for the cantilever
bondingpad  p1 [left_anchor][l=20u w=20u h=4u oz=pi]
beam3d  p1 [left_anchor middle][l=L w=w h=h  oz=angle T=T]
commonground  p1 [right_anchor][l=10u w=10u h=4u oz=0]
beam3d  p1 [middle right_anchor ][l=L w=w h=h  oz=-angle T=T]

% beam3d  p1 [middle up][l=20u w=10u h=h  oz=pi/2]
beam3d  p1 [middle dn][l=20u w=10u h=h  oz=-pi/2]

