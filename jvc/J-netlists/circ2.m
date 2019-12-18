%net=cho_load('circ2.m');[q,k]=cho_dc(net);figure(1);cho_display(net,q);
%q(lookup_coord(net,'b','z'));q50=q(1:6,1)
   
uses mumps.net
anchor p1 [a][l=10u w=10u h=10u oz=pi]
semicircularbeam p1 [a b(1)] [radius=100u w=20u h=2u alpha=pi/4 ]
K=8*3
for k=1:K-1
   [
      semicircularbeam p1 [b(k) b(k+1)] [radius=100u w=20u h=2u alpha=pi/4 oz=k*pi/4]
   ]
f3d * [b(K-1)][F=5u oy=pi/2]
f3d * [b(K-4-1)][F=5u oy=pi/2]


