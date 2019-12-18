subnet subrectangle4cornernodetags [UR UL LR LL] [ ]
[
   
%net=cho_load('mfrectangle4cornernodetags.m');figure(1);cho_display(net);
%uses mumps2.net

%temporary anchor 
%anchor p1 [a(1,1)][l=10u w=10u h=10u oz=-90] 
%upper right a(1,   wres) 
%lower right a(lres,wres) 
%lower left  a(lres,1) 
%upper left  a(1,   1) 

pi=3.141592
oy=0
oz=0
l=200u
w=200u
wresolution=5 
lresolution=5 
wres=wresolution
lres=lresolution
W=w/wres
L=l/lres

for i=1:wres [ %for each row
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
beam3d p1 [a(1,wres) UR][l=ltag w=wtag oz=theta] %upper right 
beam3d p1 [a(lres,wres) LR][l=ltag w=wtag oz=-theta] %lower right 
beam3d p1 [a(1,1) UL][l=ltag w=wtag oz=(180-theta)] %upper left  
beam3d p1 [a(lres,1) LL][l=ltag w=wtag oz=(180+theta)] %lower left 

anchor p1 [UL] [l=10u w=10u h=10u oz=-90] 
anchor p1 [LL] [l=10u w=10u h=10u oz=-90] 

f3d * [UR][ F=200u oy=-90] 
f3d * [LR][ F=200u oy=-90] 

beam3d p1 [UL B][l=l w=w oz=180]
f3d * [B][ F=400u oy=-90] 
]
