% Function used internally to display a 3D beam based on
% Hermite spline interpolants between the end points.
%
% Inputs:
%   q_global - displacements of the two nodes in global coordinates:
%              [x1; y1; z1; rx1; ry1; rz1;  x2; y2; z2; rx2; ry2; rz2]
%   R        - transform local to global coordinates
%   P       - the coordinates of the first node of the beam ("root point")
%   L, W, H  - beam length, width, and height
%
% Output:
%   Displays the beam in a Matlab plot

%command:
%q=zeros(12,1);P=[0;0;0];R=rot2global(0,0,pi/8);W=10e-6;H=2e-6;radius=100e-6;alpha=pi/4;figure(2);displaycircbeam4(q, P, R, W, H, radius, alpha);view(2);
%net=cho_load('mirror2.m');[q,k]=cho_dc(net);figure(2);cho_display(net,q);

function displaycircbeam4(q_global, P, R, W, H, radius, alpha)

%By: Jason Vaughn Clark - Oct 1998.
% Modified by David Bindel, 7/2001
% Modified by Jason Vaughn Clark, 10/2001

resolution = 80; %Plot resolution.
L=2*radius*sin(alpha/2); %x-projection of semicircle, in a bowl position.
A=[0,0,0,0;W/2,W/2,-W/2,-W/2;-H/2,H/2,H/2,-H/2]; %Beam yz-plane cross section 
s=linspace(0,L,resolution); %Sample points along beam

%Global displacements
q=q_global;
 x1G=q(1);  y1G=q(2);  z1G=q(3); 
rx1G=q(4); ry1G=q(5); rz1G=q(6); 
 x2G=q(7);  y2G=q(8);  z2G=q(9); 
rx2G=q(10);ry2G=q(11);rz2G=q(12); 

%Local displacements 
q=reshape(R'*reshape(q_global,3,4), 12,1); 
 x1L=q(1);  y1L=q(2);  z1L=q(3); 
rx1L=q(4); ry1L=q(5); rz1L=q(6); 
 x2L=q(7);  y2L=q(8);  z2L=q(9); 
rx2L=q(10);ry2L=q(11);rz2L=q(12); 

%Bowl displacements 
q=reshape(rot2global(0,0,-alpha/2)*(R'*reshape(q_global,3,4)), 12,1); 
 x1B=q(1);  y1B=q(2);  z1B=q(3); 
rx1B=q(4); ry1B=q(5); rz1B=q(6); 
 x2B=q(7);  y2B=q(8);  z2B=q(9); 
rx2B=q(10);ry2B=q(11);rz2B=q(12); 

%The warped neutral beam axis h. Add a straight beam polynomial to a semicircle.  
%hxB(s) would go from -L/2 to L/2 if xB=0.
hxB=x1B+... %Node1 position,
    (1+(x2B-x1B)/L)*s+... %to node2 position.
    -L/2; %Shift it to be symmetric about y-axis.
%A bowl shape in the xy-plane. Ends of the bowl are at y((+/-)L/2)=0. 
hyB=hermite_cubic(L,y1B,rz1B,y2B,rz2B,s)+... %Distortion of a straight beam
    sqrt(radius^2-(L/2).^2)+... %mapped onto a semicircular contour
    -sqrt(radius^2-(s-L/2).^2); %in the shape of a bowl.
%hzB(s) displacements.
hzB=hermite_cubic(L,z1B,-ry1B,z2B,-ry2B,s); %Distortion in z.
%hzB=hermite_cubic(L,z1B,ry1B,z2B,ry2B,s); %Distortion in z.

%Reset node1 of hB.
hxB=hxB-x1B; 
hyB=hyB-y1B; 
hzB=hzB-z1B; 

%String the cross sections, A, along the hB(s) beam contour. 
sigma=linspace(-alpha/2,alpha/2,resolution); %Angle from -alpha/2 to alpha/2.
S=linspace(0,radius*sin(alpha),resolution); %Sample points from 0:radius*sin(alpha).
f=asin(S./radius)/alpha; %rotation factor, 0:1.

%b1=rot2global(pi/2,0,0)*(rot2global(0,0,pi/6*0)*(rot2global(0,0,pi/6/2)*(rot2global(0,0,pi/6/2)*[0;20e-6/2;2e-6/2]))) + [100e-6*sin(pi/6);100e-6*(1-cos(pi/6));0]

for corner=1:4 %The four rectangular corners of the cross section, A.
   p=zeros(3,resolution); %Reset p, the [x;y;z] coordinates of each corner for each sample point.
   pcorner=[A(1,corner);A(2,corner);A(3,corner)]; %Choose a corner at its initial position.
   for i=1:resolution %Pick a sample point to rotate and position along the bowl contour.
      p(:,i)=rot2global(rx2G*f(i),ry2G*f(i),rz2G*f(i))*...         
             (R*... %Rotate to global.
             (rot2global(0,0,alpha/2)*... %Rotate to local.
             (rot2global(0,0,sigma(i))*... %Rotation due to bowl curvature, -alpha/2 to alpha/2.
             pcorner))); %Corner point being rotated into orientation.
      p(:,i)=rot2global(rx1G,ry1G,rz1G)*...
             (p(:,i)+... %Rotate due to displacement of node1.
             (R*... %Rotate to global   
             (rot2global(0,0,alpha/2)*... %Rotate to local.
             [hxB(i)+L/2;hyB(i);hzB(i)])));
   end
%if corner==2
%a=   p(:,1)           +P+[x1G;y1G;z1G]
%b=   p(:,resolution)  +P+[x1G;y1G;z1G]
%end

  
  
   u(corner,:) = p(1,:) + P(1) + x1G; %translate to global x positions
   v(corner,:) = p(2,:) + P(2) + y1G; %translate to global y positions
   w(corner,:) = p(3,:) + P(3) + z1G; %translate to global z positions
end

%figure(2); plot3(p1(1,:),p1(2,:),p1(3,:)); view(2);
%figure(3); plot3(p2(1,:),p2(2,:),p2(3,:)); view(2);
%figure(4); plot3(p3(1,:),p3(2,:),p3(3,:)); view(2);
%figure(5); plot3(hx+P(1),hy+P(2),hz+P(3)); view(2);
%pause

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

