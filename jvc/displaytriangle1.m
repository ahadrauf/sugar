% Function used internally to display a 3D beam based on
% Hermite spline interpolants between the end points.
%
% Inputs:
%   q_global - displacements of the two nodes in global coordinates:
%              [x1; y1; z1; rx1; ry1; rz1;  x2; y2; z2; rx2; ry2; rz2]
%   Rp       - the coordinates of the first node of the beam ("root point")
%   L, W, H  - beam length, width, and height
%   R        - rotation matrix
% Output:
%   Displays the beam in a Matlab plot

%M=10e-6;p.h=10e-6;p.L1=100e-6;p.L2=200e-6;p.ox=0;p.oy=0;p.oz=0;p.alpha=30;q=[0;0;0;0;0;0;  0;0;M;0;0;0;  0;0;0;0;0;0]; displaytriangle1(q,p); 
function displaytriangle1(q, parameter)
%By: Jason Vaughn Clark - Jul2001.

%q => [x0,y0,z0,rx0,ry0,rz0,  x1,y1,z1,rx1,ry1,rz1,  x2,y2,z2,rx2,ry2,rz2]'
%p0 => [x,y,z]'. Position of nondeflected node 1.
iresolution = 10;		%Plot resolution.
jresolution = 10;		%Plot resolution.

p1=rot2local(parameter.ox,parameter.oy,parameter.oz)*[parameter.L1;0;0];
p2=rot2local(parameter.ox,parameter.oy,parameter.oz)*[parameter.L2;0;0];
H=parameter.h;
%side1 in unrotated frame = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%input parameter
L=parameter.L1; 
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
for x = 0:L/iresolution:L
   i=i+1;
   u1(i) = ax0 + ax1.*x;
   v1(i) = ay0 + ay1.*x + ay2.*x.^2 + ay3.*x.^3;
   w1(i) = az0 + az1.*x + az2.*x.^2 + az3.*x.^3;
end

%side2 in unrotated frame  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%parameters
L=parameter.L2;
alpha=parameter.alpha*pi/180;
R=rot2local(ox,oy,oz+alpha);
%displacements
Q=R*[q(1);q(2);q(3)]; %rotate node1 xyz to local
x1=Q(1);
y1=Q(2);
z1=Q(3);
Q=R*[q(4);q(5);q(6)]; %rotate node1 oxoyoz to local
ox1=Q(1);
oy1=Q(2);
oz1=Q(3);
Q=R*[q(1+12);q(2+12);q(3+12)]; %rotate node2 xyz to local
x2=Q(1);
y2=Q(2);
z2=Q(3);
Q=R*[q(4+12);q(5+12);q(6+12)]; %rotate node2 oxoyoz to local
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
for x = 0:L/iresolution:L
   i=i+1;
   u2(i) = ax0 + ax1.*x;
   v2(i) = ay0 + ay1.*x + ay2.*x.^2 + ay3.*x.^3;
   w2(i) = az0 + az1.*x + az2.*x.^2 + az3.*x.^3;
   Q=R'*[u2(i);v2(i);w2(i)]; %rotate alpha
   u2(i) = Q(1);
   v2(i) = Q(2);
   w2(i) = Q(3);
end

%side21 in unrotated frame  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%parameters
%calc beta
A=[1;0;0]; B=[parameter.L1-parameter.L2*cos(alpha); 0-parameter.L2*sin(alpha);0];
beta=-acos(dot(A,B)/norm(A)/norm(B));
R=rot2local(ox,oy,oz+beta);   
%calc length
A=[parameter.L2*cos(alpha);parameter.L2*sin(alpha);0];
B=[parameter.L1;0;0];
L=norm(A-B);
LL=L;
%displacements
%for j=resolution:-1:1   
for j=1:jresolution
   %x,y initial points
   ix1=parameter.L1*j/jresolution;
   iy1=0;
   iz1=0;
   ix2=parameter.L2*cos(alpha)*j/jresolution;
   iy2=parameter.L2*sin(alpha)*j/jresolution;
   iz2=0;
   L=LL*j/jresolution;
   %x,y displacement points
   dx1=-(u1(j+1)-ix1);
   dy1=-(v1(j+1)-iy1);
   dz1=-(w1(j+1)-iz1);
   dx2=-(u2(j+1)-ix2);
   dy2=-(v2(j+1)-iy2);
   dz2=-(w2(j+1)-iz2);
   %line
   Q=R*[dx1;dy1;dz1];
   x1=Q(1);
   y1=Q(2);
   z1=Q(3);
   Q=j/jresolution*R*[q(4+12);q(5+12);q(6+12)]; %wieghted rotations 
   ox1=Q(1);
   oy1=Q(2);
   oz1=Q(3);
   Q=R*[dx2;dy2;dz2];
   x2=Q(1);
   y2=Q(2);
   z2=Q(3);
   Q=j/jresolution*R*[q(4+6);q(5+6);q(6+6)]; %wieghted rotations 
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
   for x = 0:L/iresolution:L
      i=i+1;
      u21(i,j) = ax0 + ax1.*x;
      v21(i,j) = ay0 + ay1.*x + ay2.*x.^2 + ay3.*x.^3; 
      w21(i,j) = az0 + az1.*x + az2.*x.^2 + az3.*x.^3; 
      Q=R'*[u21(i,j);v21(i,j);w21(i,j)];
      u21(i,j) = Q(1) + u2(j+1);
      v21(i,j) = Q(2) + v2(j+1);
      w21(i,j) = Q(3) + w2(j+1);
   end
end


%figure(1);clf;
%hold on;
%rotate3d on;
%plot3(u2,v2,w2,'b');
%plot3(u1,v1,w1,'b',u2,v2,w2,'b');
%for i=1:iresolution
%   plot3(u1,v1,w1,'b',u2,v2,w2,'b',u21(:,i),v21(:,i),w21(:,i),'b');
%end

[I,J]=size(u21);
x=zeros(I,J*2+1+2);
y=zeros(I,J*2+1+2);
z=zeros(I,J*2+1+2);

for i=1:I
   z(i,1)=H/2;
   z(i,J*2+2)=-H/2;
   z(i,J*2+3)=H/2;
end

for i=1:I
   for j=1:J
      x(i,j+1)=u21(i,j);
      y(i,j+1)=v21(i,j);
      z(i,j+1)=w21(i,j)+H/2;
   end
   for j=1:J
      x(i,1+J+j)=u21(i,J+1-j);
      y(i,1+J+j)=v21(i,J+1-j);
      z(i,1+J+j)=w21(i,J+1-j)-H/2;
   end
end

figure(1);clf;
color=pink;
color(:,2)=color(:,2)*(0.5 + 0.5*rand);
color(:,1)=color(:,1)*(0.5 + 0.5*rand);
color(:,3)=color(:,3)*(0.5 + 0.5*rand);
colormap(color);
axis vis3d;
axis equal;
%surfl([zeros(size(x,1),1),x],[zeros(size(x,1),1),y],[zeros(size(x,1),1),z]);
surfl(x,y,z);shading interp;colormap(color);shading interp;
rotate3d on;axis equal;
