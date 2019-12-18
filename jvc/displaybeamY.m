% Function used internally to display a 3D beam based on
% Hermite spline interpolants between the end points.
%
% Inputs:
%   q_global - displacements of the two nodes in global coordinates:
%              [x1; y1; z1; rx1; ry1; rz1;  x2; y2; z2; rx2; ry2; rz2]
%   Rp       - the coordinates of the first node of the beam ("root point")
%   L, W, H  - beam length, width, and height
%
% Output:
%   Displays the beam in a Matlab plot

function displaybeam(q_global, Rp, R, L, W, H)

%By: Jason Vaughn Clark - Oct 1998.

%q => [x1,y1,z1,rx1,ry1,rz1,  x2,y2,z2,rx2,ry2,rz2]'
%R => [x,y,z]'. Position of nondeflected node 1.
%parameter => Beam geometry parameters.

resolution = 20;		%Plot resolution.

%Transformation matrix from system origin to rotation R.
% TR = rot2local(Rrx,Rry,Rrz); 
TR = R';

%Transform q from global to local coordinates:
%q = rotate2local_yzx_12x12(Rrx,Rry,Rrz)*q_global;   
q = reshape(TR*reshape(q_global,3,4), 12,1);

%Starting corners of a relaxed beam, along x: 
% points A# (B#) represent corner points of node1 (node2).
A = [0    0    0    0;
     W/2  W/2 -W/2 -W/2;
    -H/2  H/2  H/2 -H/2];
B = A;
     
%Apply the displacement rotation to the beam-end corners. 
%Rotate from O to displacement, but still on x axis:
%Node 1 q-rotation.
T1 = rot2local(q(4),q(5),q(6))';  	%Node1 rotation about x,y,z.
T2 = rot2local(q(10),q(11),q(12))';	%Node2 rotation about x,y,z.
A = T1*A;
B = T2*B;

%Plot valuing: in local coords.
%x-space:
x = linspace(0,L,resolution); 

%u(x) coefficients: u(x) = ax0 + ax1.*x.
ax1 = (L+q(7)-q(1)) / L;

%v(x) coefficients. v(x) = ay0 + ay1.*x + ay2.*x.^2 + ay3.*x.^3.
ay1 = q(6);
ay2 = 3*(q(8)-q(2))/L^2	- (2*q(6)+q(12))/L;
ay3 = -2*(q(8)-q(2))/L^3 + (q(6)+q(12))/L^2;

%w(x) coefficients. w(x) = az0 + az1.*x + az2.*x.^2 + az3.*x.^3.
az1 = -q(5);
az2 = 3*(q(9)-q(3))/L^2 + (2*q(5)+q(11))/L;
az3 = -2*(q(9)-q(3))/L^3 - (q(5)+q(11))/L^2;

%Parameterized points from each corner at node1 to node2. 
% v = ay0 + ay1.*x + ay2.*x.^2 + ay3.*x.^3;
% w = az0 + az1.*x + az2.*x.^2 + az3.*x.^3;   
%A# to B# parameterized points.
for corner = 1 : 4
  
  %Calculate points
  u(corner,1:resolution) = ax1*x + ...
      (A(1,corner) + (B(1,corner) - A(1,corner))*x/L);
  v(corner,1:resolution) = ay1*x + ay2*x.^2 + ay3*x.^3 + ...
      (A(2,corner) + (B(2,corner) - A(2,corner))*x/L);
  w(corner,1:resolution) = az1*x + az2*x.^2 + az3*x.^3 + ...
      (A(3,corner) + (B(3,corner) - A(3,corner))*x/L);
  
  points = TR' * [u(corner,:); v(corner,:); w(corner,:)]; 
  u(corner,:) = points(1,:) + Rp(1) + q_global(1); 
  v(corner,:) = points(2,:) + Rp(2) + q_global(2); 
  w(corner,:) = points(3,:) + Rp(3) + q_global(3); 
 
end %for corner   

%Surfacing.

X = [u(1:4,:); u(1,:)];
Y = [v(1:4,:); v(1,:)];
Z = [w(1:4,:); w(1,:)];

surfl(X,Y,Z); 
