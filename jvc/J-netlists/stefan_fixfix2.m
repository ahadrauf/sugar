% p.E=165e9;p.w=1e-6;p.h=2e-6;I=p.w^3*p.h/12;p.F=5e-6;p.L=100e-6; net=cho_load('stefan_fixfix2.m',p); dq=cho_dc(net); Sugar_deflection=dqval(net,dq,'b','y'), ;fprintf('rotating_anchors \t(F*L^3)/(48*E*I)  = \t%g\n',(p.F*p.L^3)/(48*p.E*I));fprintf('fixfixed_anchors \t(F*L^3)/(192*E*I) = \t%g\n',(p.F*p.L^3)/(192*p.E*I));figure(1);cho_display(net,dq)

param E=165e9
param w=1u
param h=2u
param I=w^3*h/12
param F=5u
param L=100u
param angle=-90

uses mumps.net
uses stdlib.net

anchorrot  p1 [a]   [l=5u w=4u h=h oz=deg(angle)]
beam2d     p1 [a b] [l=L/2 w=w h=h Youngsmodulus=E]
beam2d     p1 [b c] [l=L/2 w=w h=h Youngsmodulus=E]
anchorrot  p1 [c]   [l=5u w=4u h=h oz=deg(angle)]
f3d        *  [b]   [F=F oz=deg(angle)]

