% Display the mechanical structure described by a netlist.
%
% Inputs:
%   net - the netlist structure
%   q   - (Optional) - the displacement of the original structure
%         If q is unspecified, the undisplaced structure will
%         be displayed.

function cho_display(net,q,options)

clf;
grid on;
hold on;
view(0,90);
xlabel('X - horizontal  [m]');		%x-axis label.
ylabel('Y - vertical  [m]');		%y-axis label.
zlabel('Z - out of plane  [m]');	%z-axis label.

if nargin == 1
    q = [];
    disp=0; %PM - to identify disp or no disp case
else
    disp=1; %PM - to identify disp or no disp case
end

%green (in 30s)
%yellow (in 40s)
%blue (in 1s)
%red in 100s
color0=35;
color1=0;  
color2=90;
color3=40;

for i = 1 : length(net.elements)
    
    elt = net.elements(i);
    
    if isfield(elt.parameter,'layer')
        switch elt.parameter.layer
            case 'p0'
                color=color0;
            case 'p1'
                color=color1;
            case 'p2'
                color=color2;
            case 'met'
                color=color3;
            otherwise
                color=color1; %default
        end
    else
        continue %double check this
    end
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
    %  eltmodel
    feval( elt.model, 'display', elt.R, elt.parameter, qlocal, 0, net.nodes(elt.node_ids), color, disp );
    
end

rotate3d on;
shading interp;
axis equal;
axis vis3d;

if disp==1
    colorbar
    title('Total displacement [m]');
end

if nargin==3
    if isfield(options,'shownodes')
        if strcmp(options.shownodes,'on')
            shownodes(net);
        end
    end
end

hold off;
