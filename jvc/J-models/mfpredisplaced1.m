subnet predisplacedbeam [a(1) a(11)][l=* w=* qx=* qy=* qz=* qox=* qoy=* qoz=*]
resolution=10
x1=0
y1=0
z1=0
ox1=0
oy1=0
oz1=0
x2=qx
y2=qy
z2=qz
ox2=qox
oy2=qoy
oz2=qoz
L=l
%slope deflection coefficients
ax0=0
ax1=(L+x2-x1)/L
ay0=0
ay1=oz1
ay2=3*(y2-y1)/L^2-(2*oz1+oz2)/L
ay3=-2*(y2-y1)/L^3+(oz1+oz2)/L^2
az0=0
az1=-oy1
az2=3*(z2-z1)/L^2+(2*oy1+oy2)/L
az3=-2*(z2-z1)/L^3-(oy1+oy2)/L^2
%segmented beams
for i=1:resolution
   x=L*i/resolution 
   dvdx=ay1+2*ay2.*x+3*ay3.*x.^2 
   dwdx=az1+2*az2.*x+3*az3.*x.^2 
   beam3d parent [a(i) a(i+1)][l=L/resolution w=W ox=ox2*i/resolution oy=dvdx oz=dwdx] 
end   
