% Display the mechanical structure described by a netlist.
%
% Inputs:
%   net - the netlist structure
%   q   - (Optional) - the displacement of the original structure
%         If q is unspecified, the undisplaced structure will
%         be displayed.

function cho_display(net,q,options);

clf;
%grid on;
hold on;
%view(0,90);
xlabel('X - horizontal  [m]');		%x-axis label.
ylabel('Y - vertical  [m]');		%y-axis label.
zlabel('Z - out of plane  [m]');	%z-axis label.

if nargin == 1
  q = [];
end

for i = 1 : length(net.elements)

  elt = net.elements(i);

  if isempty(elt.var_ids)
    qlocal = [];
  else
    j = find(elt.var_ids ~= 0);
    jdx = elt.var_ids(j);
    qlocal = zeros(length(elt.var_ids),1);
    if (~isempty(q))
      qlocal(j) = q(jdx);
    end
  end
  
  feval( elt.model, 'display', elt.R, elt.parameter, qlocal, 0, ...
         net.nodes(elt.node_ids) );

end

%rotate3d on;
color=pink;
color(:,2)=color(:,2)*(0.5 + 0.5*rand);
color(:,1)=color(:,1)*(0.5 + 0.5*rand);
color(:,3)=color(:,3)*(0.5 + 0.5*rand);
colormap(color);
shading interp;
%axis equal;
%axis vis3d;

if nargin==3
   if isfield(options,'shownodes')
      if strcmp(options.shownodes,'on')
         shownodes(net);
      end
   end
end

hold off;
