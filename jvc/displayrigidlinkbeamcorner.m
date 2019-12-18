% Function used internally to display a 3D beam with cloaked links.
%
% Inputs:
%   q_global - displacements of the two link nodes in global coordinates:
%              [x1; y1; z1; rx1; ry1; rz1;  x2; y2; z2; rx2; ry2; rz2]
%   R2glo    - transform local to global coordinates
%   Rp       - the coordinates of the first node of the beam ("root point")
%   param    - l,w,h,L1,L2,etc
%
% Output:
%   Displays the shifted beam in a Matlab plot

function displayrigidlinkbeamcorner(q_linkglo, Rp, R2glo, param)

%By: Jason Vaughn Clark - Oct 1998.
% Modified by David Bindel, 7/2001
% Modified by Jason Clark, 10/2001

%parameters
L=param.l; %beam parameters
W=param.w;
H=param.h;
resolution=20; %Plot resolution.
if isfield(param,'L1'),L1=param.L1;else,L1=0;end
if isfield(param,'L2'),L2=param.L2;else,L2=0;end
s=linspace(-2*L1,L+L2*2,resolution); 
if ~isfield(param,'oy1'),oy1=0;else, oy1=param.oy1; end
if ~isfield(param,'oy2'),oy2=0;else, oy2=param.oy2; end
if ~isfield(param,'oz1'),oz1=0;else, oz1=param.oz1; end
if ~isfield(param,'oz2'),oz2=0;else, oz2=param.oz2; end
a1=rot2local(0,oy1,oz1)'*[L1;0;0];  
Tau1=[0,-a1(3),a1(2);a1(3),0,-a1(1);-a1(2),a1(1),0];
a2=-rot2local(0,oy2,oz2)'*[L2;0;0];
Tau2=[0,-a2(3),a2(2);a2(3),0,-a2(1);-a2(2),a2(1),0];
i3=eye(3);
o3=zeros(3);
Tau1=[i3,o3;Tau1,i3];
Tau2=[i3,o3;Tau2,i3];
o6=zeros(6);
TAU=[Tau1,o6;o6,Tau2];

%local beam displacements
q_linkloc=reshape(R2glo'*reshape(q_linkglo,3,4), 12,1); %Transform global link to local link
q_beamloc=TAU'*q_linkloc; %Transform local link to local beam
 x1 = q_beamloc(1);   y1 = q_beamloc(2);   z1 = q_beamloc(3);
rx1 = q_beamloc(4);  ry1 = q_beamloc(5);  rz1 = q_beamloc(6);
 x2 = q_beamloc(7);   y2 = q_beamloc(8);   z2 = q_beamloc(9);
rx2 = q_beamloc(10); ry2 = q_beamloc(11); rz2 = q_beamloc(12);

%Beam cross seciton
A = [0    0    0    0;
     W/2  W/2 -W/2 -W/2;
    -H/2  H/2  H/2 -H/2];
B = A;
A = rot2local(rx1, ry1, rz1)' * A;
B = rot2local(rx2, ry2, rz2)' * B;

% -- Build the coefficients for Hermite interpolants through
%    corresponding corners.  The result will be a displacement
%    field (hx, hy, hz) to be added to the undisplaced beam shape.
%    Just stretch along x axis, use cubic interpolants for y and z
%    Now a cubic interpolants for y displacement
hx = x1 + (1 + (x2-x1)/L)*s;
hy = hermite_cubic(L, y1,rz1, y2,rz2, s);
hz = hermite_cubic(L, z1,-ry1, z2,-ry2, s);

% -- Compute interpolants through the corners and put them into the global coordinate frame
Rb=Rp+R2glo*rot2local(0,oy1,oz1)'*[L1;0;0]; 
for corner = 1 : 4
  px = A(1,corner) + (B(1,corner)-A(1,corner))*s/L + hx;
  py = A(2,corner) + (B(2,corner)-A(2,corner))*s/L + hy;
  pz = A(3,corner) + (B(3,corner)-A(3,corner))*s/L + hz;
  p = R2glo * [px; py; pz];
  u(corner,:) = p(1,:) + Rb(1);
  v(corner,:) = p(2,:) + Rb(2);
  w(corner,:) = p(3,:) + Rb(3);
end

% -- Plot the surfaces
X = [u(1:4,:); u(1,:)];
Y = [v(1:4,:); v(1,:)];
Z = [w(1:4,:); w(1,:)];
surfl(X,Y,Z); 

% Evaluate a cubic Hermite interpolant with data given at 0 and L.
% Use Newton's divided difference form (see your favorite intro
%  numerical analysis book).
function [f] = hermite_cubic(L, f0, f00, fL, fLL, s)
f0L   = (fL - f0)/L;
f00L  = (f0L - f00)/L;
f0LL  = (fLL - f0L)/L;
f00LL = (f0LL - f00L)/L;
f = f0 + s.*(f00 + s.*(f00L + (s-L).*f00LL));
