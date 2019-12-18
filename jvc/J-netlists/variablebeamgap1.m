%

i=0;
for V=4:2:10
    %update the parameter
        p.V=V; 
        net=cho_load('variablebeamgap1.m',p);
    %display results
        q=cho_dc(net); 
        i=i+1; figure(i); 
        cho_display(net,q); 
    %display deflection meters
    y=q(lookup_coord(net,'c','y'));
end

uses mumps.net
uses stdlib.net

param V=10 %Voltage parameter input. If there is no input parameter, then it defaults to V=10.
Vsrc     *  [A f] [V=V] %The parameter is used here.

%everything else is the same.
eground  *  [f]   []
anchor   p1 [A]   [l=5u w=10u  oz=deg(180)]
beam2de  p1 [A b] [l=100u w=2u h=2u oz=0 R=100]
gap2de   p1 [b c D E] [l=100u w1=10u w2=5u oz=0 gap=2u R1=100 R2=100]
eground  *  [D]   []
anchor   p1 [D]   [l=5u w=10u oz=-deg(90)]
anchor   p1 [E]   [l=5u w=10u oz=-deg(90)]
eground  *  [E]   []

