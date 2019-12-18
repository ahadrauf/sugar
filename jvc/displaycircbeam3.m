% Function used internally to display a 3D beam based on
% Hermite spline interpolants between the end points.
%
% Inputs:
%   q_global - displacements of the two nodes in global coordinates:
%              [x1; y1; z1; rx1; ry1; rz1;  x2; y2; z2; rx2; ry2; rz2]
%   R        - transform local to global coordinates
%   Rp       - the coordinates of the first node of the beam ("root point")
%   L, W, H  - beam length, width, and height
%
% Output:
%   Displays the beam in a Matlab plot

function displaycircbeam3(q_global, Rp, R, W, H, radius, alpha)

%By: Jason Vaughn Clark - Oct 1998.
% Modified by David Bindel, 7/2001
% Modified by Jason Vaughn Clark, 10/2001

resolution = 80; %Plot resolution.
L=radius*sin(alpha); %x-projection of semicircle
R2loc1=R'; %rotate to local coordinates due to node 1
R2glo1=R; %rotate to global coordinates due to node 1
R2alpha=rot2local(0,0,alpha)'; %rotate to alpha

%Transform q to local coordinates 
q = q_global;
 x1g = q(1);  y1g =  q(2);  z1g =  q(3); 
 rx1g = q(4); ry1g =  q(5); rz1g =  q(6); 
 x2g = q(7);  y2g =  q(8);  z2g =  q(9); 
rx2g = q(10);ry2g = q(11); rz2g = q(12); 

q = reshape(R2loc1*reshape(q_global,3,4), 12,1); 
 x1 = q(1);  y1 =  q(2);  z1 =  q(3); 
 rx1 = q(4); ry1 =  q(5); rz1 =  q(6); 
 x2 = q(7);  y2 =  q(8);  z2 =  q(9); 
rx2 = q(10);ry2 = q(11); rz2 = q(12); 

%rotations
R2glo2=rot2local(rx2, ry2, rz2)'; %rotate to global coordinates due to node 2

%Compute the displacement of the corner points A & B for node 1 and 2. Then apply rotational displacements.
A = [0,0,0,0;W/2,W/2,-W/2,-W/2;-H/2,H/2,H/2,-H/2]; 

%Where along the beam will we put sample points / knots?
s = linspace(0,L,resolution); 

%The warped neutral beam axis. 
%Add a straight beam polynomial to a semicircle. Step along x axis. Cubic interpolants for y(x) and z(x). 
hx = x1 + (1 + (x2-x1)/L)*s; 
hy = hermite_cubic(L,y1,rz1,y2,rz2,s) + radius-sqrt(radius^2-s.^2); 
hz = hermite_cubic(L,z1,-ry1,z2,-ry2,s) +0*hermite_cubic(L,0,rx1,0,rx2,radius-sqrt(radius^2-s.^2));

%String A-to-B cross sections along the beam.
f=asin(s./radius)/alpha;
for corner = 1 : 4 %The four rectangular corners
  for i=1:resolution
     px = A(1,corner); %x coords of a corner
     py = A(2,corner); %y coords of a corner
     pz = A(3,corner); %z coords of a corner     
     R2alpha=rot2local(0,0,alpha*f(i))'; %rotate to alpha
     R2glo1=rot2local(rx1g*(1-f(i)),ry1g*(1-f(i)),rz1g*(1-f(i)))'; %rotate to global coordinates due to node 1
     R2glo2=rot2local(rx2g*f(i),ry2g*f(i),rz2g*f(i))'; %rotate to global coordinates due to node 2     
     p(:,i)=R2glo1*R2glo2*R*R2alpha*[px;py;pz] + R*[hx(i);hy(i);hz(i)];
  end  
  u(corner,:) = p(1,:) + Rp(1); %translate to spacial x positions
  v(corner,:) = p(2,:) + Rp(2); %translate to spacial y positions
  w(corner,:) = p(3,:) + Rp(3); %translate to spacial z positions
end

%i=20
%R 
%R2alpha=rot2local(0,0,alpha*f(i))' 
%R2glo1=rot2local(rx1*(1-f(i)),ry1*(1-f(i)),rz1*(1-f(i)))' 
%R2glo2=rot2local(rx2*f(i),ry2*f(i),rz2*f(i))' 

%i=1
R2alpha=rot2local(0,0,alpha*f(i))'; %rotate to alpha
     R2glo1=rot2local(rx1*(1-f(i)),ry1*(1-f(i)),rz1*(1-f(i)))'; %rotate to global coordinates due to node 1
     R2glo2=rot2local(rx2*f(i),ry2*f(i),rz2*f(i))'; %rotate to global coordinates due to node 2     
%     R2glo2,R2glo1,R,R2alpha
%r=     R2glo2*R2glo1*R*R2alpha
%i=20     
     R2alpha=rot2local(0,0,alpha*f(i))'; %rotate to alpha
     R2glo1=rot2local(rx1*(1-f(i)),ry1*(1-f(i)),rz1*(1-f(i)))'; %rotate to global coordinates due to node 1
     R2glo2=rot2local(rx2*f(i),ry2*f(i),rz2*f(i))'; %rotate to global coordinates due to node 2     
%     R2glo2,R2glo1,R,R2alpha
%     r=R2glo2*R2glo1*R*R2alpha
%rx1,ry1,rz1     
%rx2,ry2,rz2     
     
%Plot the surfaces
X = [u(1:4,:); u(1,:)];
Y = [v(1:4,:); v(1,:)];
Z = [w(1:4,:); w(1,:)];
surfl(X,Y,Z); 

%Evaluate a cubic Hermite interpolant with data given at 0 and L. Use Newton's divided difference form.
function [f] = hermite_cubic(L, f0, f00, fL, fLL, s)
f0L   = (fL - f0)/L;
f00L  = (f0L - f00)/L;
f0LL  = (fLL - f0L)/L;
f00LL = (f0LL - f00L)/L;
f = f0 + s.*(f00 + s.*(f00L + (s-L).*f00LL));

