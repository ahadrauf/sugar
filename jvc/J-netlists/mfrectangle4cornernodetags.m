%mfrectangle4cornernode.m
%net=cho_load('mfrectangle4cornernodetags.m');q=cho_dc(net);figure(1);cho_display(net,q);
%q(lookup_coord(net,'UR','y')),q(lookup_coord(net,'B','y'))  

uses mumps2.net

anchor p1 [UL] [l=10u w=10u h=10u oz=-90] 
anchor p1 [LL] [l=10u w=10u h=10u oz=-90] 
f=1000u

%temporary anchor 
%anchor p1 [a(1,1)][l=10u w=10u h=10u oz=-90] 
%upper right a(1,   wres) 
%lower right a(lres,wres) 
%lower left  a(lres,1) 
%upper left  a(1,   1) 

pi=3.141592
param ox=0, oy=0, oz=0
param l=400u, w=200u, wresolution=20 lresolution=20 
wres=wresolution
lres=lresolution
W=w/wres
L=l/lres

for i=1:wres [ %for each row
   anchor p1 [a(i,1) cc(i)][l=10u w=10u h=10u oz=-90]
   for j=1:lres-1 [ %for each row segment
      beam3d p1 [a(i,j) a(i,j+1)][l=L w=W] 
   ] 
] 
for i=1:lres [ %for each column
   for j=1:wres-1 [ %for each column segment
      beam3d p1 [a(j,i) a(j+1,i)][l=W w=L oz=-90] 
   ] 
] 

ltag=sqrt((W*W/4 + L*L/4))
wtag=ltag
theta=atan(W/L)*180/pi 
ht=2u
beam3d p1 [a(1,wres) UR][l=ltag w=wtag oz=theta h=ht] %upper right 
beam3d p1 [a(lres,wres) LR][l=ltag w=wtag oz=-theta h=ht] %lower right 
beam3d p1 [a(1,1) UL][l=ltag w=wtag oz=(180-theta) h=ht] %upper left  
beam3d p1 [a(lres,1) LL][l=ltag w=wtag oz=(180+theta) h=ht] %lower left 

f3d * [UR][ F=f/2 oy=-90] 
f3d * [LR][ F=-f/2 oy=-90] 

beam3d p1 [UL BB][l=w/2 w=2u oz=-90] 
anchor p1 [BB][l=10u w=10u h=10u oz=-90] 
%beam3d p1 [BB B][l=l-(l/lres) w=w oz=180] 
beam3d p1 [BB B][l=l w=w oz=180] 
%f3d * [B][ F=f oz=-90] 
f3d * [B][ M=-f/2*w/2 oz=0] 


%beam3d p1 [UL c1][l=l/4 w=w oz=180]
%beam3d p1 [c1 c2][l=l/4 w=w oz=180]
% beam3d p1 [c1 c12][l=l/4 w=w oz=90]
%beam3d p1 [c2 c3][l=l/4 w=w oz=180]
% beam3d p1 [c3 c13][l=l/4 w=w oz=90]
%beam3d p1 [c3 B][l=l/4 w=w oz=180]

