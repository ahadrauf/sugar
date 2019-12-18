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

function displaycircbeam2(q_global, Rp, R, L, W, H, radius, alpha)

%By: Jason Vaughn Clark - Oct 1998.
% Modified by David Bindel, 7/2001
% Modified by Jason Vaughn Clark - Sep 2001

resolution = 10;		%Plot resolution.

% -- Transform q to local coordinates

q = reshape(R'*reshape(q_global,3,4), 12,1);

 x1 =  q(1);  y1 =  q(2);  z1 =  q(3);
rx1 =  q(4); ry1 =  q(5); rz1 =  q(6);
 x2 =  q(7);  y2 =  q(8);  z2 =  q(9);
rx2 = q(10); ry2 = q(11); rz2 = q(12);


% -- Compute the displacement of the corner points
%    (A for the corners at one end, B for the corners at the other)

A = [0    0    0    0;
     W/2  W/2 -W/2 -W/2;
    -H/2  H/2  H/2 -H/2];
B = A;


% -- Apply the (small) deflections to the corners due to
%    rotational displacement
     
A = rot2local(rx1, ry1, rz1)' * A;
B = rot2local(rx2, ry2, rz2)' * B;


% -- Where along the beam will we put sample points / knots?

s = linspace(0, L, resolution); 


% -- Build the coefficients for Hermite interpolants through
%    corresponding corners.  The result will be a displacement
%    field (hx, hy, hz) to be added to the undisplaced beam shape.
%    Just stretch along x axis, use cubic interpolants for y and z
%    Now a cubic interpolants for y displacement

hx = x1 + (1 + (x2-x1)/L)*s;
hy = circularshape(radius,alpha,y1,rz1,y2,rz2,s);
hz = hermite_cubic(L, z1,-ry1, z2,-ry2, s);

% -- Compute interpolants through the corners and put them into the
%    global coordinate frame

for corner = 1 : 4
   
for i=1:resolution
  px(i) = A(1,corner) + (B(1,corner)-A(1,corner))*s(i)/L ;
  py(i) = A(2,corner) + (B(2,corner)-A(2,corner))*s(i)/L ;
  pz(i) = A(3,corner) + (B(3,corner)-A(3,corner))*s(i)/L ;  
  p = rot2local(0,0,alpha*s(i)/L)'*[px(i); py(i); pz(i)];
  px(i)=p(1)+ hx(i);
  py(i)=p(2)+ hy(i);
  pz(i)=p(3)+ hz(i);
end
  
  p = R * [px; py; pz];

  u(corner,:) = p(1,:) + Rp(1);
  v(corner,:) = p(2,:) + Rp(2);
  w(corner,:) = p(3,:) + Rp(3);
 
end

% -- Plot the surfaces

X = [u(1:4,:); u(1,:)];
Y = [v(1:4,:); v(1,:)];
Z = [w(1:4,:); w(1,:)];

surfl(X,Y,Z); 


% ---------
% Evaluate a cubic Hermite interpolant with data given at 0 and L.
% Use Newton's divided difference form (see your favorite intro
%  numerical analysis book).

function [f] = hermite_cubic(L, f0, f00, fL, fLL, s)
f0L   = (fL - f0)/L;
f00L  = (f0L - f00)/L;
f0LL  = (fLL - f0L)/L;
f00LL = (f0LL - f00L)/L;
f = f0 + s.*(f00 + s.*(f00L + (s-L).*f00LL));

function [cb]=circularshape(radius,alpha,y1,oz1,y2,oz2,x)
L=radius*sin(alpha);
q1=y1;
q2=tan(oz1);
%q2=tan(q2); %angle to slope
q3=y2;
q4=tan(oz2);
%q4=tan(q4); %angle to slope
psi1=1-3*(x./L).^2+2*(x./L).^3;
psi2=x.*(1-x./L).^2;
psi3=3*(x./L).^2-2*(x./L).^3;
psi4=x.^2./L.*(x./L-1);
Y = psi1.*q1 + psi2.*q2 + psi3.*q3 + psi4.*q4;
cb=(radius-sqrt(radius.^2-x.^2))+Y;
axis equal

