%mfrectangle4cornernode.m
%net=cho_load('mfrectangle4cornernode.m');figure(1);cho_display(net);
uses mumps2.net

%temporary anchor 
anchor p1 [a(1,1)][l=10u w=10u h=10u oz=-90] 
%upper right a(1,     wres+1) 
%lower right a(lres+1,wres+1) 
%lower left  a(lres+1,1) 
%upper left  a(1,     1) 

pi=3.141592
param ox=0, oy=0, oz=0
param l=200u, w=200u, wresolution=20 lresolution=20 
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

