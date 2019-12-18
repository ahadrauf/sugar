
subnet perforated_subnet2 [leftA leftB rightA rightB rightC][numberholes=* w=* h=* l=* whorizontal=* wvertical=*] 
[
   %LH=(l-wvertical)/numberholes-wvertical
   LH=14u
   %LV=w-whorizontal*2
%   LV=14u
   LV=10u
   %H=0.98
   H=1
   %Make n-1 inner holes
   for i=1:numberholes-1
   [
      beam3d parent [a(i) a(i+1)][l=LH w=whorizontal h=h ] %Horizontal.
      beam3d parent [b(i) b(i+1)][l=LH w=whorizontal h=h ] %Horizontal.
      beam3dlink parent [a(i+1) b(i+1)][l=LV w=wvertical h=h*H L1=whorizontal/2 L2=whorizontal/2 oz=-pi/2] %Vertical.
%      beam3dlink parent [a(i+1) b(i+1)][l=LV w=wvertical h=h oz=-pi/2] %Vertical.
   ]

   %Make left end
      beam3dlinkcorner parent [leftA  a(1)][l=LH w=whorizontal h=h ] %Horizontal.
      beam3dlinkcorner parent [leftB  b(1)][l=LH w=whorizontal h=h ] %Horizontal.
      beam3dlink parent [leftA  leftB][l=LV w=wvertical h=h*H L1=whorizontal/2 L2=whorizontal/2 oz=-pi/2] %Vertical.
      beam3dlink parent [a(1)  b(1)][l=LV w=wvertical h=h*H L1=whorizontal/2 L2=whorizontal/2 oz=-pi/2] %Vertical extra one.
%      beam3dlink parent [leftA  leftB][l=LV w=wvertical h=h*H oz=-pi/2] %Vertical.
%      beam3dlink parent [a(1)  b(1)][l=LV w=wvertical h=h*H oz=-pi/2] %Vertical extra one.

   %beam3dlinkcorner parent [a(1) left][l=LV/2 w=wvertical h=h L1=whorizontal/2 oz=-pi/2] %Vertical.
   %beam3dlinkcorner parent [left b(1)][l=LV/2 w=wvertical h=h L2=whorizontal/2 oz=-pi/2] %Vertical.
   
   %Make right end
   beam3dlinkcorner parent [a(numberholes) rightA][l=LH w=whorizontal h=h ] %Horizontal.
   beam3dlinkcorner parent [b(numberholes) rightB][l=LH w=whorizontal h=h ] %Horizontal.
   beam3dlink parent [rightA rightC][l=LV/2 w=wvertical h=h*H L1=whorizontal/2 oz=-pi/2] %Vertical.
   beam3dlink parent [rightC rightB][l=LV/2 w=wvertical h=h*H L2=whorizontal/2 oz=-pi/2] %Vertical.
%   beam3dlink parent [rightA rightC][l=LV/2 w=wvertical h=h*H oz=-pi/2] %Vertical.
%   beam3dlink parent [rightC rightB][l=LV/2 w=wvertical h=h*H oz=-pi/2] %Vertical.
]


