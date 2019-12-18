%net=cho_load('perforated4.m');figure(1);cho_display(net);
%net=cho_load('perforated4.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q); 
%net=cho_load('perforated4.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q); a=abs(q(130)),c=sqrt(2*(9e-6)^2);x=c*cos(pi/4+a);y=c*sin(pi/4+a);xy=9e-6;d=sqrt((xy-x)^2 + (xy-y)^2)
%net=cho_load('perforated4.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q); a=abs(q(length(q)-2)),c=sqrt(2*(9e-6)^2);x=c*cos(pi/4+a);y=c*sin(pi/4+a);xy=9e-6;d=sqrt((xy-x)^2 + (xy-y)^2)

uses mumps3.net
uses perforated_subnet2.m

L=1u

anchor p1 [leftA][l=L w=L h=L oz=-pi]
anchor p1 [leftB][l=7u w=7u h=7u oz=pi/2]

perforated_subnet2 p1 [leftA leftB rightA rightB rightC][numberholes=9 w=18u h=18u l=140u whorizontal=4u wvertical=4u] 

beam3d p1 [rightA xa][l=(9u*9)-2u w=4u h=18u oz=pi/2 ]
beam3d p1 [rightB xb][l=(9u*9)-2u w=4u h=18u oz=-pi/2 ]

F = (1e-1) * (1/4) * (1/2) /(1)/1e2
%F = (1e-1) * (1/4) * (1/2) * (9/7)/1e3/2

f3d * [xa][F=F oy=pi/2]
f3d * [xb][F=-F oy=pi/2]
%f3d * [rightC][M=2*F*9e-6]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%beam3d p1 [leftA x][l=100u w=10u h=2u ]
%f3d * [x][M=2*F*2.5e-6 / 10]

