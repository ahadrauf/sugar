% net=cho_load('net_ucb_contact_zip_1.m');[q]=cho_dc(net);figure(1);cho_display(net,q); q(lookup_coord(net,'B','z'))
% net=cho_load('net_ucb_contact_zip_1.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('net_ucb_contact_zip_1.m');figure(1);cho_display(net);


uses mumps.net
param V=0
 
anchor      p1 [e]      [l=5u w=10u  oz=deg(180)]
beam2de_b   p1 [e a(1)]    [l=50u w=2u h=2u oz=0 Rs=10]

for i = 1 : 4
    [
        gap2de_h    p1 [a(i) a(i+1) b(i) b(i+1)][l=50u/4 w1=2u w2=5u oz=0 gap=4u Rs=10 R_insulator=100]
        anchor      p1 [b(i)]      [l=5u w=10u oz=-deg(90)]
        anchor      p1 [b(i+1)]    [l=5u w=10u oz=-deg(90)]
        eground     *  [b(i)]   []
        eground     *  [b(i+1)] []
	]

Vsrc     *  [e g] [V=V]
eground  *  [g]   []


