%[ux,vx,wx,rxx,ryx,rzx]=beamdeflection(net,q,name)
%net is from net=cho_load('netlist');
%q is from q=cho_dc(net);
%name is from the name given to the element in the 'netlist' such as 
%    anchor     p1 [A]   [l=10u w=10u]
%    hi  beam3d p1 [A B] [l=100u w=6u h=2u oz=0]
%    bye beam3d p1 [B C] [l=100u w=6u h=2u oz=45]
%    f3d        *  [C] [F=100u oz=-90]
%    where the first beam is named "hi" and the second beam element is named "bye."
%From the command line this whole netlist run would be like the following
%   net=cho_load('netlist');q=cho_dc(net);figure(1);cho_display(net,q);
%   beamdeflection(net,q,'hi');
%You may need to increase the size of the figure window for clarity.

function [ux,vx,wx,rxx,ryx,rzx]=beamdeflection(net,q,name)

%find the index of the element name
for i=1:length(net.elements)
   if strcmp(net.elements(i).name,name)
      elementindex=i;
      break;
   end
end

%find node names
node1=net.elements(elementindex).node{1};
node2=net.elements(elementindex).node{2};

%set up element geometric quantities
q_global=zeros(12,1);
if lookup_coord(net,node1,'x')
   q_global(1)=q(lookup_coord(net,node1,'x'));
   q_global(2)=q(lookup_coord(net,node1,'y'));
   q_global(3)=q(lookup_coord(net,node1,'z'));
   q_global(4)=q(lookup_coord(net,node1,'rx'));
   q_global(5)=q(lookup_coord(net,node1,'ry'));
   q_global(6)=q(lookup_coord(net,node1,'rz'));
end
if lookup_coord(net,node2,'x')
   q_global(7)=q(lookup_coord(net,node2,'x'));
   q_global(8)=q(lookup_coord(net,node2,'y'));
   q_global(9)=q(lookup_coord(net,node2,'z'));
   q_global(10)=q(lookup_coord(net,node2,'rx'));
   q_global(11)=q(lookup_coord(net,node2,'ry'));
   q_global(12)=q(lookup_coord(net,node2,'rz'));
end
Rp=[0;0;0]; %postion of node 1 is set to the origin.
R=net.elements(elementindex).R; %rotation matrix
L=net.elements(elementindex).parameter.l; %length
W=net.elements(elementindex).parameter.w; %width
H=net.elements(elementindex).parameter.h; %layer thickness
resolution = 100;		%Plot resolution.
TR = R';
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

xx=0:L/resolution:L;
ux=q(1)+ax1.*xx;
vx=q(2)+ay1.*xx + ay2.*xx.^2 + ay3.*xx.^3;
wx=q(3)+az1.*xx + az2.*xx.^2 + az3.*xx.^3;
rxx=q(4)+(L+q(10)-q(4))/L.*xx;
ryx=q(5)+(L+q(11)-q(5))/L.*xx;
rzx=q(6)+(L+q(12)-q(6))/L.*xx;
figure;clf;
subplot(7,1,1);plot(xx,ux);ylabel('x');title('Beam displacesments x,y,z,ox,oy,oz as funcitons of axial lenght l=0:L.');grid on;
subplot(7,1,2);plot(xx,vx);ylabel('y');grid on;
subplot(7,1,3);plot(xx,wx);ylabel('z');grid on;
subplot(7,1,4);plot(xx,rxx);ylabel('ox');grid on;
subplot(7,1,5);plot(xx,ryx);ylabel('oy');grid on;
subplot(7,1,6);plot(xx,rzx);ylabel('oz');grid on;
subplot(7,1,7);
grid on;
hold on;
view(0,90);
rotate3d on;
axis equal;
xlabel('x          beam picture'); %x-axis label.
ylabel('y'); %y-axis label.
zlabel('z'); %z-axis label.

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

rotate3d on;
color=pink;
color(:,2)=color(:,2)*(0.5 + 0.5*rand);
color(:,1)=color(:,1)*(0.5 + 0.5*rand);
color(:,3)=color(:,3)*(0.5 + 0.5*rand);
colormap(color);
shading interp;
axis equal;
axis vis3d;
hold off;
