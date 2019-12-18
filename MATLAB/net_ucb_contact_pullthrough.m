% net=cho_load('net_ucb_contact_pullthrough.m');[q]=cho_dc(net);figure(1);cho_display(net,q); q(lookup_coord(net,'B','z'))
% net=cho_load('net_ucb_contact_pullthrough.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('net_ucb_contact_pullthrough.m');figure(1);cho_display(net);


uses mumps.net
param V=0
L1 = 2u
h = 2u
anchor      p1 [e]      [l=10u w=10u  oz=deg(180)]
beam2de_b   p1 [e b1]   [l=50u w=2u h=2u oz=0 Rs=10]
beam2de_b   p1 [b1 b2]  [l=L1 w=2u h=h oz=deg(-90) Rs=10]
%beam2de_b   p1 [b1 b22]  [l=L1*2 w=2u h=h oz=deg(0) Rs=10]
%beam2de_b   p1 [b2 b3]  [l=L1 w=2u h=2u oz=deg(0) Rs=10]
beam2de_b   p1 [b3 a(1)]  [l=L1 w=2u h=h oz=deg(90) Rs=10]

gap2de_h_passive    p1 [b2 b3 bb2 bb3][l=L1 w1=2u w2=5u oz=0 gap=4u-L1 Rs=10 R_insulator=100 h=h]
gap2de_h_passive    p1 [c3 c2 cc2 cc3][l=L1 w1=2u w2=5u oz=0 gap=4u-L1 Rs=10 R_insulator=100 h=h]
anchor      p1 [bb3] [l=5u w=10u oz=-deg(90)]
anchor      p1 [bb2] [l=5u w=10u oz=-deg(90)]
eground     *  [bb2] []
eground     *  [bb3] []
anchor      p1 [cc3] [l=5u w=10u oz=-deg(90)]
anchor      p1 [cc2] [l=5u w=10u oz=-deg(90)]
eground     *  [cc2] []
eground     *  [cc3] []

N = 6
for i = 1 : N
    [
        gap2de_h    p1 [a(i) a(i+1) b(i) b(i+1)][l=70u/N w1=2u w2=5u oz=0 gap=4u Rs=10 R_insulator=100 h=h]
        anchor      p1 [b(i)]      [l=5u w=10u oz=-deg(90)]
        anchor      p1 [b(i+1)]    [l=5u w=10u oz=-deg(90)]
        eground     *  [b(i)]   []
        eground     *  [b(i+1)] []
	]

beam2de_b   p1 [a(N+1) c3]  [l=L1 w=2u h=h oz=deg(-90) Rs=10]
%beam2de_b   p1 [c3 c2]  [l=L1 w=2u h=2u oz=deg(0) Rs=10]
beam2de_b   p1 [c2 c1]  [l=L1 w=2u h=h oz=deg(90) Rs=10]
beam2de_b   p1 [c1 c0]  [l=L1 w=2u h=h oz=deg(0) Rs=10]

Vsrc     *  [e g] [V=V]
eground  *  [g]   []


