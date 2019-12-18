%subnet_comb for 5 or more fingers.
%This subnet makes 1 half of a comb drive. 
%Jason Vaughn Clark - Nov2001 

subnet subnet_comb [left middle right][numberfingers=* gap=* fingerwidth=* fingerheight=* fingerlength=* supportwidth=* supportheight=* V=*] 
[
   fingerrigidlink=supportwidth/2
   supportrigidlink=fingerwidth/2
   fingerspacing=gap*2+fingerwidth %
   %Make the left (n-5)/2 inner fingers
   for i=2:(numberfingers-5)/2+1
   [
      beam3dlinkcorner parent [a(i) a(i+1)][l=fingerspacing w=supportwidth h=supportheight L1=fingerwidth/2 L2=fingerwidth/2] %Horizontal support; growing from left to right.
      beam3dlink parent [a(i) b(i)][l=fingerlength w=fingerwidth h=fingerheight oz=-pi/2 L1=supportwidth/2] %Fingers hang down from support.
   ]
   %Make the right (n-5)/2 inner fingers
   for i=(numberfingers-5)/2+4:numberfingers-2
   [
      beam3dlinkcorner parent [a(i) a(i+1)][l=fingerspacing w=supportwidth h=supportheight L1=fingerwidth/2 L2=fingerwidth/2] %Horizontal support; growing from left to right.
      beam3dlinkcorner parent [a(i+1) b(i)][l=fingerlength w=fingerwidth h=fingerheight oz=-pi/2 L1=supportwidth/2] %Fingers hang down from support.
   ]
   %Make the left side  
   beam3dlinkcorner parent [left a(2)][l=fingerspacing w=supportwidth h=supportheight L1=fingerwidth/2 L2=fingerwidth/2] %Horizontal support; growing from left to right.
   beam3dlink parent [left bleft][l=fingerlength w=fingerwidth h=fingerheight oz=-pi/2 L1=supportwidth/2] %Fingers hang down from support.
   %Make the middle
   beam3dlinkcorner parent [a((numberfingers-5)/2+2) middle][l=fingerspacing w=supportwidth h=supportheight L1=fingerwidth/2 L2=fingerwidth/2] %Horizontal support; growing from left to right.
   beam3dlink parent [a((numberfingers-5)/2+2) bml][l=fingerlength w=fingerwidth h=fingerheight oz=-pi/2 L1=supportwidth/2] %Middle left finger.   
   beam3dlink parent [middle bmiddle][l=fingerlength w=fingerwidth h=fingerheight oz=-pi/2 L1=supportwidth/2] %Middle finger
   beam3dlinkcorner parent [middle a((numberfingers-5)/2+4)][l=fingerspacing w=supportwidth h=supportheight L1=fingerwidth/2 L2=fingerwidth/2] %Horizontal support; growing from left to right.
   beam3dlink parent [a((numberfingers-5)/2+4) bmr][l=fingerlength w=fingerwidth h=fingerheight oz=-pi/2 L1=supportwidth/2] %Middle right finger.   
   %Make the right side  
   beam3dlinkcorner parent [a(numberfingers-1) right][l=fingerspacing w=supportwidth h=supportheight L1=fingerwidth/2 L2=fingerwidth/2] %Horizontal support; growing from left to right.
   beam3dlink parent [right bright][l=fingerlength w=fingerwidth h=fingerheight oz=-pi/2 L1=supportwidth/2] %Fingers hang down from support.
   %Forces and moments   
   fingerforce=(8.854e-12)*2*fingerheight*V^2/gap/2 %Force on a rotor finger due to two adjacent stator fingers.
   L=(fingerspacing+fingerwidth)*(numberfingers-1)-fingerwidth/2 %Length of the support, from the middle of the left finger to the middle of the right finger.
   P=fingerforce*numberfingers/L %Newtons/length. 
   F1=P*L/2 %Force on node left.
   F6=(P*L^2)/12 %Moment on node left.
   F7=P*L/2 %Force on node right.
   F12=-(P*L^2)/12 %Moment on node right.
   f3d * [left][F=F1 oz=-pi/2] %Force on node left.
   f3d * [left][M=F6 oy=-pi/2] %Moment on node left.
   f3d * [right][F=F7 oz=-pi/2] %Force on node right.
   f3d * [right][M=F12 oy=-pi/2] %Moment on node right.
]


