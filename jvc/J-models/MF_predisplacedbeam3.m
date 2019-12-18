%net=cho_load('predisp1.m');
%figure(1);cho_display(net);
function [output] = MF_predisplacedbeam3(flag, R, params, q, t, nodes, varargin);
switch(flag)
case 'vars'
   output.dynamic={1{'x','y','z','rx','ry','rz'};2{'x','y','z','rx','ry','rz'}}; 
case 'M' 
   global recycle_m_predisplaced
   num_m=length(recycle_recycle_m_predisplacedbeam2); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_m%Go through each one
      if (recycle_m_predisplacedbeam2{i}.w==params.w & ... 
         recycle_m_predisplacedbeam2{i}.qox1==params.qox1 & ... 
         recycle_m_predisplacedbeam2{i}.qoy1==params.qoy1 & ... 
         recycle_m_predisplacedbeam2{i}.qoz1==params.qoz1 & ... 
         recycle_m_predisplacedbeam2{i}.qox2==params.qox2 & ... 
         recycle_m_predisplacedbeam2{i}.qoy2==params.qoy2 & ... 
         recycle_m_predisplacedbeam2{i}.qoz2==params.qoz2 & ... 
         recycle_m_predisplacedbeam2{i}.qx2==params.qx2 & ... 
         recycle_m_predisplacedbeam2{i}.qy2==params.qy2 & ... 
         recycle_m_predisplacedbeam2{i}.qz2==params.qz2 & ... 
         recycle_m_predisplacedbeam2{i}.l==params.l & ... 
         recycle_m_predisplacedbeam2{i}.h==params.h) 
         found_it=i; %Store the cell index of the right one found.
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated k.
      m=recycle_m_predisplacedbeam2{found_it}.m; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('predisplacedbeam2.m',params); %Load the equivalent discretized model. 
      m=assemble_system(net,'M',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      m=reordermatrixtop(m,node1ids,node2ids); %Reorder the matrix st.
      [m]=matrixcondensationtop(m); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_m_predisplacedbeam2{found_it+1}.qox1==params.qox1;
      recycle_m_predisplacedbeam2{found_it+1}.qoy1==params.qoy1;
      recycle_m_predisplacedbeam2{found_it+1}.qoz1==params.qoz1;
      recycle_m_predisplacedbeam2{found_it+1}.qox2==params.qox2;
      recycle_m_predisplacedbeam2{found_it+1}.qoy2==params.qoy2;
      recycle_m_predisplacedbeam2{found_it+1}.qoz2==params.qoz2;
      recycle_m_predisplacedbeam2{found_it+1}.qx2==params.qx2;
      recycle_m_predisplacedbeam2{found_it+1}.qy2==params.qy2;
      recycle_m_predisplacedbeam2{found_it+1}.qz2==params.qz2;
      recycle_m_predisplacedbeam2{found_it+1}.l==params.l;
      recycle_m_predisplacedbeam2{found_it+1}.h==params.h;
   end
   %Rotate m into the global frame.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   RL=rigidlinkmatrix(params); %Compute the rigid link transformations.
   output=R*(RL*m*RL')*R'; %Return the global 12*12 stiffness matrix.
case 'D'
   global recycle_d_predisplaced
   num_d=length(recycle_recycle_d_predisplacedbeam2); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_d%Go through each one
      if (recycle_d_predisplacedbeam2{i}.w==params.w & ... 
         recycle_d_predisplacedbeam2{i}.qox1==params.qox1 & ... 
         recycle_d_predisplacedbeam2{i}.qoy1==params.qoy1 & ... 
         recycle_d_predisplacedbeam2{i}.qoz1==params.qoz1 & ... 
         recycle_d_predisplacedbeam2{i}.qox2==params.qox2 & ... 
         recycle_d_predisplacedbeam2{i}.qoy2==params.qoy2 & ... 
         recycle_d_predisplacedbeam2{i}.qoz2==params.qoz2 & ... 
         recycle_d_predisplacedbeam2{i}.qx2==params.qx2 & ... 
         recycle_d_predisplacedbeam2{i}.qy2==params.qy2 & ... 
         recycle_d_predisplacedbeam2{i}.qz2==params.qz2 & ... 
         recycle_d_predisplacedbeam2{i}.l==params.l & ... 
         recycle_d_predisplacedbeam2{i}.h==params.h) 
         found_it=i; %Store the cell index of the right one found.
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated k.
      d=recycle_d_predisplacedbeam2{found_it}.d; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('predisplacedbeam2.m',params); %Load the equivalent discretized model. 
      d=assemble_system(net,'D',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      d=reordermatrixtop(d,node1ids,node2ids); %Reorder the matrix st.
      [d]=matrixcondensationtop(d); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_d_predisplacedbeam2{found_it+1}.qox1==params.qox1;
      recycle_d_predisplacedbeam2{found_it+1}.qoy1==params.qoy1;
      recycle_d_predisplacedbeam2{found_it+1}.qoz1==params.qoz1;
      recycle_d_predisplacedbeam2{found_it+1}.qox2==params.qox2;
      recycle_d_predisplacedbeam2{found_it+1}.qoy2==params.qoy2;
      recycle_d_predisplacedbeam2{found_it+1}.qoz2==params.qoz2;
      recycle_d_predisplacedbeam2{found_it+1}.qx2==params.qx2;
      recycle_d_predisplacedbeam2{found_it+1}.qy2==params.qy2;
      recycle_d_predisplacedbeam2{found_it+1}.qz2==params.qz2;
      recycle_d_predisplacedbeam2{found_it+1}.l==params.l;
      recycle_d_predisplacedbeam2{found_it+1}.h==params.h;
   end
   %Rotate d into the global frame.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   RL=rigidlinkmatrix(params); %Compute the rigid link transformations.
   output=R*(RL*d*RL')*R'; %Return the global 12*12 stiffness matrix.
case 'K'
   global recycle_k_predisplaced
   num_k=length(recycle_recycle_k_predisplacedbeam2); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_k%Go through each one
      if (recycle_k_predisplacedbeam2{i}.w==params.w & ... 
         recycle_k_predisplacedbeam2{i}.qox1==params.qox1 & ... 
         recycle_k_predisplacedbeam2{i}.qoy1==params.qoy1 & ... 
         recycle_k_predisplacedbeam2{i}.qoz1==params.qoz1 & ... 
         recycle_k_predisplacedbeam2{i}.qox2==params.qox2 & ... 
         recycle_k_predisplacedbeam2{i}.qoy2==params.qoy2 & ... 
         recycle_k_predisplacedbeam2{i}.qoz2==params.qoz2 & ... 
         recycle_k_predisplacedbeam2{i}.qx2==params.qx2 & ... 
         recycle_k_predisplacedbeam2{i}.qy2==params.qy2 & ... 
         recycle_k_predisplacedbeam2{i}.qz2==params.qz2 & ... 
         recycle_k_predisplacedbeam2{i}.l==params.l & ... 
         recycle_k_predisplacedbeam2{i}.h==params.h) 
         found_it=i; %Store the cell index of the right one found.
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated k.
      k=recycle_k_predisplacedbeam2{found_it}.k; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('predisplacedbeam2.m',params); %Load the equivalent discretized model. 
      k=assemble_system(net,'K',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      k=reordermatrixtop(k,node1ids,node2ids); %Reorder the matrix st.
      [k]=matrixcondensationtop(k); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_k_predisplacedbeam2{found_it+1}.qox1==params.qox1;
      recycle_k_predisplacedbeam2{found_it+1}.qoy1==params.qoy1;
      recycle_k_predisplacedbeam2{found_it+1}.qoz1==params.qoz1;
      recycle_k_predisplacedbeam2{found_it+1}.qox2==params.qox2;
      recycle_k_predisplacedbeam2{found_it+1}.qoy2==params.qoy2;
      recycle_k_predisplacedbeam2{found_it+1}.qoz2==params.qoz2;
      recycle_k_predisplacedbeam2{found_it+1}.qx2==params.qx2;
      recycle_k_predisplacedbeam2{found_it+1}.qy2==params.qy2;
      recycle_k_predisplacedbeam2{found_it+1}.qz2==params.qz2;
      recycle_k_predisplacedbeam2{found_it+1}.l==params.l;
      recycle_k_predisplacedbeam2{found_it+1}.h==params.h;
   end
   %Rotate k into the global frame.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   RL=rigidlinkmatrix(params); %Compute the rigid link transformations.
   output=R*(RL*k*RL')*R'; %Return the global 12*12 stiffness matrix.
case 'pos'
   params
   pause
   net=cho_load('predisplacedbeam2.m',params); %Load the equivalent discretized model. 
   net
   pause
%   s=length(net.nodes);
%   pos=net.nodes(s).pos;
   % Compute relative positions of beam nodes
%   [RL1,RL2]=rigidlinkcoordinates(params,R); %[x;y;z] coordinated of both rigidlinks.
%   P=RL2-RL1; %difference of rigidlink coordinates.
pos=[0;0;0];
P=pos;
   output=R*[0,pos(1)+P(1);0,pos(2)+P(2);0,pos(3)+P(3)]; %R*[node1postions,node2postions]   
case 'display'
   L=params.l+params.qx2-params.qx1;
   [RL1,RL2]=rigidlinkcoordinates(params,R); %[x;y;z] coordinated of both rigidlinks.
   q(1:3,1)=q(1:3,1)+RL1; %Translation of node1 due to link1.
   q(7:9,1)=q(7:9,1)+RL1; %Translation of node2 due to link1.
   q(4)=q(4)+params.qox1;
   q(5)=q(5)+params.qoy1;
   q(6)=q(6)+params.qoz1;
   q(7)=q(7)+params.qx2;
   q(8)=q(8)+params.qy2;
   q(9)=q(9)+params.qz2;
   q(10)=q(10)+params.qox2;
   q(11)=q(11)+params.qoy2;
   q(12)=q(12)+params.qoz2;
   displaybeamcorner(q,nodes(1).pos,R,params);
case 'check'
   S=[];
   if ~isfield(params,'h'),S=strcat(S,' predisplacedbeam2 model needslayer thickness (h).');end
   if ~isfield(params,'w'),S=strcat(S,' predisplacedbeam2 model needs width (w).');end
   if ~isfield(params,'l'),S=strcat(S,' predisplacedbeam2 model length (l).');end
   if ~isfield(params,'Youngsmodulus'),S=strcat(S,' predisplacedbeam2 model needsprocess file Youngs modulus (Youngsmodulus).');end
   if ~isfield(params,'fluid'),S=strcat(S,' predisplacedbeam2 model needsprocess file viscous fluid layer thickness (fluid).');end
   if ~isfield(params,'viscosity'),S=strcat(S,' predisplacedbeam2 model needsprocess file viscosity of fluid layer (viscosity).');end
   if ~isfield(params,'density'),S=strcat(S,' predisplacedbeam2 model needsprocess file material density (density).');end
   output=[S];
end

