% net=cho_load('net_ucb_contact_1a.m');[q]=cho_dc(net);figure(1);cho_display(net,q); q(lookup_coord(net,'B','z'))
% net=cho_load('net_ucb_contact_1a.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('net_ucb_contact_1a.m');figure(1);cho_display(net);


uses mumps.net
param V=0
 
anchor      p1 [e]      [l=5u w=10u  oz=deg(180)]
beam2de_b   p1 [e a]    [l=100u w=2u h=2u oz=0 Rs=10]
gap2de_f    p1 [a b c d][l=100u w1=10u w2=5u oz=0 gap=2u Rs=10 dielectric=10]
anchor      p1 [c]      [l=5u w=10u oz=-deg(90)]
anchor      p1 [d]      [l=5u w=10u oz=-deg(90)]


Vsrc     *  [e g] [V=V]
eground  *  [g]   []
eground  *  [c]   []
eground  *  [d]   []

%f2d * [d] [F=100000n oz=-pi/2]

