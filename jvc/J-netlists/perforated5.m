%net=cho_load('perforated4.m');figure(1);cho_display(net);
%net=cho_load('perforated4.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q); 
%net=cho_load('perforated4.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q); a=abs(q(130)),c=sqrt(2*(9e-6)^2);x=c*cos(pi/4+a);y=c*sin(pi/4+a);xy=9e-6;d=sqrt((xy-x)^2 + (xy-y)^2)
%net=cho_load('perforated5.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q); 

uses mumps3.net
uses perforated_subnet6.m

L=1u

perforated_subnet6 p1 [ aa bb A B ][numberholes=10 w=18u h=18u l=140u whorizontal=4u wvertical=4u] 

%beam3dlink p1 [A xa][l=(9u*9)+0u w=4u h=18u oz=pi/2 L1=2u]
%beam3dlink p1 [B xb][l=(9u*9)+0u w=4u h=18u oz=-pi/2 L1=2u]

F = 0*(1e-1) * (1/4) * (1/2) /(1)/1e2
F = (1e-1) * (1/4) * (1/2) /(1)/1e2 *(9/7)

f3d * [A][F=F oy=pi/2]
f3d * [B][F=-F oy=pi/2]

anchor p1 [aa][l=1u w=1u h=1u oz=-pi]
anchor p1 [bb][l=1u w=1u h=1u oz=-pi]
%anchor p1 [A][l=1u w=1u h=1u oz=-pi]
%anchor p1 [B][l=1u w=1u h=1u oz=-pi]
