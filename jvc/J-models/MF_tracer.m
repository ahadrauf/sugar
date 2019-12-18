function [output] = MF_tracer(flag, R, params, q, t, nodes, varargin);

switch(flag)
case 'vars'
  output.dynamic = {1 {'x' 'y' 'z' 'rx' 'ry' 'rz'};2 {'x' 'y' 'z' 'rx' 'ry' 'rz'}};

case 'pos'
   % Compute relative positions of beam nodes
   x=params.l; %Relative x position of node2 before rotation.
   [RL1,RL2]=rigidlinkcoordinates(params); %[x;y;z] coordinated of both rigidlinks.
   P=RL2-RL1; %difference of rigidlink coordinates.
   output=R*[0,x+P(1);0,P(2);0,P(3)]; %R*[node1postions,node2postions]
case 'display'
   [RL1,RL2,r1,r2]=rigidlinkcoordinates(params); %[x;y;z] coordinated of both rigidlinks.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations. To.
   q=To'*q; %displacement of unlinked nodes.
   Rp=(nodes(1).pos-R*RL1);     
   displaybeam(q,Rp,R,params.l,params.w,params.h,varargin{1},varargin{2} ); %PM (varargin{1}and{2})
case 'check'
   S=[];
   if ~isfield(params,'h'),S=strcat(S,' beam3dlink model needslayer thickness (h).');end
   if ~isfield(params,'w'),S=strcat(S,' beam3dlink model needs width (w).');end
   if ~isfield(params,'l'),S=strcat(S,' beam3dlink model needsangle sweep (alpha).');end
   if ~isfield(params,'Youngsmodulus'),S=strcat(S,' beam3dlink model needsprocess file Youngs modulus (Youngsmodulus).');end
   if ~isfield(params,'fluid'),S=strcat(S,' beam3dlink model needsprocess file viscous fluid layer thickness (fluid).');end
   if ~isfield(params,'viscosity'),S=strcat(S,' beam3dlink model needsprocess file viscosity of fluid layer (viscosity).');end
   if ~isfield(params,'density'),S=strcat(S,' beam3dlink model needsprocess file material density (density).');end
   output=[S];
case 'writebin'
  fid = varargin{1};
  var_ids = varargin{2};
  if (fid < 0)
    output = 1;
  else
    writebeam(fid, nodes(1).pos, R, params.l, params.w, params.h, var_ids);
  end
otherwise
  output = [];
end

