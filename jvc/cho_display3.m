% Display the mechanical structure described by a netlist.
%
% Inputs:
%   net - the netlist structure
%   q   - (Optional) - the displacement of the original structure
%         If q is unspecified, the undisplaced structure will
%         be displayed.

function cho_display3

clf;
grid on;
hold on;
view(0,90);
xlabel('X - horizontal  [m]');		%x-axis label.
ylabel('Y - vertical  [m]');		%y-axis label.
zlabel('Z - out of plane  [m]');	%z-axis label.

  x=2e-5;
q1=[0;0;0;  0;0;0];
q2=[0;0;0;  0;0;0]; 
q3=[0;0;x;  0;0;0];

  q=[q1; q2];
  P=[0;0;0];
  R=rot2global(0,0,0);
  W=20e-6;H=2e-6;radius=100e-6;alpha=pi/6;
  displaycircbeam4(q, P, R, W, H, radius, alpha);
  
  q=[q2 ;q3];
  P=[radius*sin(alpha);radius*(1-cos(alpha));0];
  R=rot2global(0,0,pi/6);
  W=20e-6;H=2e-6;radius=100e-6;alpha=pi/6;
  displaycircbeam4(q, P, R, W, H, radius, alpha);
  
  
      
rotate3d on;
color=pink;
color(:,2)=color(:,2)*(0.5 + 0.5*rand);
color(:,1)=color(:,1)*(0.5 + 0.5*rand);
color(:,3)=color(:,3)*(0.5 + 0.5*rand);
colormap(color);
shading interp;
axis equal;
axis vis3d;

if nargin==3
   if isfield(options,'shownodes')
      if strcmp(options.shownodes,'on')
         shownodes(net);
      end
   end
end

hold off;
