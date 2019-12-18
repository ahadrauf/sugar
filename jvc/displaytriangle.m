% Inputs:
%   q        - displacements of the 3 nodes in global coordinates.
%              [x1; y1; z1; rx1; ry1; rz1;     x2; y2; z2; rx2; ry2; rz2;     x3; y3; z3; rx3; ry3; rz3]
%   pos      - the coordinates of the first node of the beam ("root point")
%   L1,L2,H  - beam length, width, and height
%   R        - rotation matrix
% Output:
%   Displays the beam in a Matlab plot

function displaytriangle(q,pos,parameter,R)
%By: Jason Vaughn Clark - Jul2001.

%q => [x0,y0,z0,rx0,ry0,rz0,  x1,y1,z1,rx1,ry1,rz1,  x2,y2,z2,rx2,ry2,rz2]'
%p0 => [x,y,z]'. Position of nondeflected node 1.
resolution = 20;		%Plot resolution.

%side1 in unrotated frame = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%input parameter
L=parameter.l1; 
ox=parameter.ox;oy=parameter.oy;oz=parameter.oz; 
%displacements
R=rot2local(ox,oy,oz); 
Q=R*[q(1);q(2);q(3)]; %rotate node1 displacements to local
x1=Q(1);
y1=Q(2);
z1=Q(3);
Q=R*[q(4);q(5);q(6)]; 
ox1=Q(1);
oy1=Q(2);
oz1=Q(3);
Q=R*[q(7);q(8);q(9)]; %rotate node2 local
x2=Q(1);
y2=Q(2);
z2=Q(3);
x12_2=x2;y12_2=y2; %save this for angle
Q=R*[q(10);q(11);q(12)]; %rotate node1 oxoyoz to local
ox2=Q(1);
oy2=Q(2);
oz2=Q(3);
%slope deflection coefficients
ax0=0;
ax1=(L+x2-x1) / L;
ay0=0;
ay1=oz1;
ay2=3*(y2-y1)/L^2	- (2*oz1+oz2)/L;
ay3=-2*(y2-y1)/L^3 + (oz1+oz2)/L^2;
az0=0;
az1=-oy1;
az2=3*(z2-z1)/L^2 + (2*oy1+oy2)/L;
az3=-2*(z2-z1)/L^3 - (oy1+oy2)/L^2;
%points along the side
i=0;
for x = 0:L/resolution:L
   i=i+1;
   u1(i) = ax0 + ax1.*x;
   v1(i) = ay0 + ay1.*x + ay2.*x.^2 + ay3.*x.^3;
   w1(i) = az0 + az1.*x + az2.*x.^2 + az3.*x.^3;
end

%side2 in unrotated frame  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%parameters
L=parameter.l2;
angle=parameter.angle*pi/180;
R=rot2local(ox,oy,oz+angle);
%displacements
Q=R*[q(1);q(2);q(3)]; %rotate node1 xyz to local
x1=Q(1);
y1=Q(2);
z1=Q(3);
Q=R*[q(4);q(5);q(6)]; %rotate node1 oxoyoz to local
ox1=Q(1);
oy1=Q(2);
oz1=Q(3);
Q=R*[q(7+6);q(8+6);q(9+6)]; %rotate node2 xyz to local
x2=Q(1);
y2=Q(2);
z2=Q(3);
x12_1=x2;y12_1=y2; %save this for angle
Q=R*[q(10+6);q(11+6);q(12+6)]; %rotate node2 oxoyoz to local
ox2=Q(1);
oy2=Q(2);
oz2=Q(3);
%slope deflection coefficients
ax0=0;
ax1=(L+x2-x1) / L;
ay0=0;
ay1=oz1;
ay2=3*(y2-y1)/L^2	- (2*oz1+oz2)/L;
ay3=-2*(y2-y1)/L^3 + (oz1+oz2)/L^2;
az0=0;
az1=-oy1;
az2=3*(z2-z1)/L^2 + (2*oy1+oy2)/L;
az3=-2*(z2-z1)/L^3 - (oy1+oy2)/L^2;
%coords
i=0;
for x = 0:L/resolution:L
   i=i+1;
   u2(i) = ax0 + ax1.*x;
   v2(i) = ay0 + ay1.*x + ay2.*x.^2 + ay3.*x.^3;
   w2(i) = az0 + az1.*x + az2.*x.^2 + az3.*x.^3;
   Q=rot2local(0,0,angle)'*[u2(i);v2(i);w2(i)]; %rotate angle
   u2(i) = Q(1);
   v2(i) = Q(2);
   w2(i) = Q(3);
end

%side12 in unrotated frame  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%parameters
angle=-atan((y12_2-y12_1)/(x12_2-x12_1)); %oz of the beam12 
R=rot2local(ox,oy,oz+angle);   
L=sqrt((u1(resolution)-u2(resolution))^2+(v1(resolution)-v2(resolution))^2+(w1(resolution)-w2(resolution))^2);
%displacements
for j=resolution:-1:1   
   x1=u2(j); %top node 
   y1=v2(j);
   z1=w2(j);
   Q=j/resolution*R*[q(10);q(11);q(12)]; %wieghted rotations 
   ox1=Q(1);
   oy1=Q(2);
   oz1=Q(3);
   x2=u1(j); %bottom node 
   y2=v1(j);
   z2=w1(j);
   Q=j/resolution*R*[q(10+6);q(11+6);q(12+6)]; %wieghted rotations 
   ox2=Q(1);
   oy2=Q(2);
   oz2=Q(3);
   %slope deflection coefficients
   ax0=0;
   ax1=(L+x2-x1) / L;
   ay0=0;
   ay1=oz1;
   ay2=3*(y2-y1)/L^2	- (2*oz1+oz2)/L;
   ay3=-2*(y2-y1)/L^3 + (oz1+oz2)/L^2;
   az0=0;
   az1=-oy1;
   az2=3*(z2-z1)/L^2 + (2*oy1+oy2)/L;
   az3=-2*(z2-z1)/L^3 - (oy1+oy2)/L^2;
   %coords
   i=0;
   for x = 0:L/resolution:L
      i=i+1;
      u12(i,j) = ax0 + ax1.*x;
      v12(i,j) = ay0 + ay1.*x + ay2.*x.^2 + ay3.*x.^3;
      w12(i,j) = az0 + az1.*x + az2.*x.^2 + az3.*x.^3;
   end
end


