% net=cho_load('torsion3.m');[q]=cho_dc(net);figure(1);cho_display(net,q);
% net=cho_load('torsion3.m');[q]=cho_dc(net);figure(1);cho_display(net,q);qy=q(lookup_coord(net,'B','y')),qz=q(lookup_coord(net,'B','z'))
% net=cho_load('torsion3.m');figure(1);cho_display(net);

%m=0;i=0;MMax=0.035e-6;for M=0:MMax/10:MMax, i=i+1;m(i)=M;p.m=M;net=cho_load('torsion3.m',p);[q]=cho_dc(net);figure(2);cho_display(net,q);qy(i)=q(lookup_coord(net,'B1','y'));,end
%figure(3);plot(qy,ox);xlabel('Lateral displacement (meters)');ylabel('Tip angle (radians)');title('constant lateral-force; increasing axial-moment'); 

f=6000u
param m = 0.035u

uses mumpsx.net
anchor p1 [A]   [l=20u w=20u oz=180 h=4u]
beam3d p1 [A B2] [l=100u w=15u h=2u]
%torsionalbeam p1 [A B1] [l=100u w=15u h=2u]

%f3d * [B1][F=6000u oz=90] 
%f3d * [B1][M=m oz=0] 
f3d * [B2][F=6000u oz=90] 
f3d * [B2][M=m oz=0] 


