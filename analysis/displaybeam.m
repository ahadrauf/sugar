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


function displaybeam(q_global, Rp, R, L, W, H,color,disp)

%By: Jason Vaughn Clark - Oct 1998.
% Modified by David Bindel, 7/2001
%  Changes by Prabhakar, March 10th 2011

resolution = 20;		%Plot resolution.


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
hy = hermite_cubic(L, y1,rz1, y2,rz2, s);
hz = hermite_cubic(L, z1,-ry1, z2,-ry2, s);


% when disp is zero PM
if disp==1
    x11=0;y11=0;z11=0;rx11=0;ry11=0;rz11=0;
    x22=0;y22=0;z22=0;rx22=0;ry22=0;rz22=0;
    hx1 = x11 + (1 + (x22-x11)/L)*s;
    hy1 = hermite_cubic(L, y11,rz11, y22,rz22, s);
    hz1 = hermite_cubic(L, z11,-ry11, z22,-ry22, s);
end  

% -- Compute interpolants through the corners and put them into the
%    global coordinate frame

for corner = 1 : 4

  px = A(1,corner) + (B(1,corner)-A(1,corner))*s/L + hx;
  py = A(2,corner) + (B(2,corner)-A(2,corner))*s/L + hy;
  pz = A(3,corner) + (B(3,corner)-A(3,corner))*s/L + hz;
  
  p = R * [px; py; pz];

  u(corner,:) = p(1,:) + Rp(1);
  v(corner,:) = p(2,:) + Rp(2);
  w(corner,:) = p(3,:) + Rp(3);
  
  if disp==1 %PM
      px1 = A(1,corner) + (B(1,corner)-A(1,corner))*s/L + hx1;
      py1 = A(2,corner) + (B(2,corner)-A(2,corner))*s/L + hy1;
      pz1 = A(3,corner) + (B(3,corner)-A(3,corner))*s/L + hz1;
      
      p1 = R * [px1; py1; pz1];
      
      u1(corner,:) = p1(1,:) + Rp(1);
      v1(corner,:) = p1(2,:) + Rp(2);
      w1(corner,:) = p1(3,:) + Rp(3);
  end
 
end

% -- Plot the surfaces

X = [u(1:4,:); u(1,:)];
Y = [v(1:4,:); v(1,:)];
Z = [w(1:4,:); w(1,:)];

h=surf(X,Y,Z);

%Diff. b/w current coordinates and coordinates when disp==0, gives the
%CData index PM
if disp==1 %Colormap indexed with displacement
    X1 = [u1(1:4,:); u1(1,:)];
    Y1 = [v1(1:4,:); v1(1,:)];
    Z1 = [w1(1:4,:); w1(1,:)];
    CData=(Y-Y1)+(X-X1)+(Z-Z1); %Total displacement
    set(h,'CData',CData);
else %Different layers with diff. colors
    set(h,'CDataMapping','direct');
    CData=zeros(size(Z));
    CData(:)=color;
    set(h,'CData',CData);
end

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
