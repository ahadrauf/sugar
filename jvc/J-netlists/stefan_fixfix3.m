%net=cho_load('stefan_fixfix3.m'); dq=cho_dc(net); Sugar_deflection=dqval(net,dq,'c','y'), figure(1);cho_display(net,dq);e=8.854e-12,w=1e-6,h=2e-6,L=100e-6,V=10,A=L*h,gap=1.5e-6,F_elec=(e*A*V^2)/(2*gap^2), E=165e9,I=w^3*h/12,fprintf('Y_analytical=(F*L^3)/(384*E*I)  = \t%g\n',(F_elec*L^3)/(384*E*I));

uses mumps.net
uses stdlib.net

param Vin = 10
L=100u
w=1u
h=1u
angle=-90
gap=1.5u
Vsrc     *  [c f] [V=Vin] 
eground  *  [f]   []
anchor   p1 [a]   [l=5u w=4u oz=deg(angle) h=h]  
gap2de   p1 [a c K Y] [l=L/2 w1=w w2=w h=h oz=0 gap=gap R1=1 R2=1 ]
gap2de   p1 [c b Y D] [l=L/2 w1=w w2=w h=h oz=0 gap=gap R1=1 R2=1 ]
anchor   p1 [b]   [l=5u w=4u oz=deg(angle) h=h]    
eground  *  [K]   []
anchor   p1 [K]   [l=5u w=10u oz=deg(angle) h=h]   
anchor   p1 [D]   [l=5u w=10u oz=deg(angle) h=h]
eground  *  [D]   []
anchor   p1 [Y]   [l=5u w=10u oz=deg(angle) h=h]
eground  *  [Y]   []
