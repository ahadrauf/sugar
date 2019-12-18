function [output] = MF_semicircularbeam(flag, R, params, q, t, nodes, varargin);
switch(flag)
case 'vars'
   output.dynamic={ 1 {'x','y','z','rx','ry','rz'}; 2 {'x','y','z','rx','ry','rz'} }; 
case 'M' 
   global recycle_m_semicircularbeamdiscrete
   num_m=length(recycle_m_semicircularbeamdiscrete); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_m%Go through each one
      if (recycle_m_semicircularbeamdiscrete{i}.w==params.w & ... %Check width.
         recycle_m_semicircularbeamdiscrete{i}.alpha==params.alpha & ... %Check angle sweep.
         recycle_m_semicircularbeamdiscrete{i}.radius==params.radius & ... %Check radius.
         recycle_m_semicircularbeamdiscrete{i}.h==params.h) %Check layer thickness.
         found_it=i; %Store the cell index of the right one found.
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated k.
      m=recycle_m_semicircularbeamdiscrete{found_it}.m; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('semicircularbeamdiscrete.m',params); %Load the equivalent discretized model. 
      m=assemble_system(net,'M',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      m=reordermatrixtop(m,node1ids,node2ids); %Reorder the matrix st.
      [m]=matrixcondensationtop(m); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_m_semicircularbeamdiscrete{found_it+1}.m=m; %Store m.
      recycle_m_semicircularbeamdiscrete{found_it+1}.w=params.w; %Store the w parameter for that m.
      recycle_m_semicircularbeamdiscrete{found_it+1}.h=params.h; %Store the h parameter for that m.
      recycle_m_semicircularbeamdiscrete{found_it+1}.radius=params.radius; %Store the radius parameter for that m.
      recycle_m_semicircularbeamdiscrete{found_it+1}.alpha=params.alpha; %Store the alpha parameter for that m.
   end
   %Rotate m into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*m*R')*To'; %Return the global 12*12 stiffness matrix.
case 'D'
   global recycle_d_semicircularbeamdiscrete
   num_d=length(recycle_d_semicircularbeamdiscrete); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_d%Go through each one
      if (recycle_d_semicircularbeamdiscrete{i}.w==params.w & ... %Check width.
         recycle_d_semicircularbeamdiscrete{i}.alpha==params.alpha & ... %Check angle sweep.
         recycle_d_semicircularbeamdiscrete{i}.radius==params.radius & ... %Check radius.
         recycle_d_semicircularbeamdiscrete{i}.h==params.h) %Check layer thickness.
         found_it=i; %Store the cell index of the right one found.
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated d.
      d=recycle_d_semicircularbeamdiscrete{found_it}.d; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('semicircularbeamdiscrete.m',params); %Load the equivalent discretized model. 
      d=assemble_system(net,'D',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      d=reordermatrixtop(d,node1ids,node2ids); %Reorder the matrix st.
      [d]=matrixcondensationtop(d); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_d_semicircularbeamdiscrete{found_it+1}.d=d; %Store d.
      recycle_d_semicircularbeamdiscrete{found_it+1}.w=params.w; %Store the w parameter for that d.
      recycle_d_semicircularbeamdiscrete{found_it+1}.h=params.h; %Store the h parameter for that d.
      recycle_d_semicircularbeamdiscrete{found_it+1}.radius=params.radius; %Store the radius parameter for that d.
      recycle_d_semicircularbeamdiscrete{found_it+1}.alpha=params.alpha; %Store the alpha parameter for that d.
   end
   %Rotate d into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*d*R')*To'; %Return the global 12*12 stiffness matrix.
case 'K'
   global recycle_k_semicircularbeamdiscrete
   num_k=length(recycle_k_semicircularbeamdiscrete); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_k%Go through each one
      if (recycle_k_semicircularbeamdiscrete{i}.w==params.w & ... %Check width.
         recycle_k_semicircularbeamdiscrete{i}.alpha==params.alpha & ... %Check angle sweep.
         recycle_k_semicircularbeamdiscrete{i}.radius==params.radius & ... %Check radius.
         recycle_k_semicircularbeamdiscrete{i}.h==params.h) %Check layer thickness.
         found_it=i; %Store the cell index of the right one found.
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated k.
      k=recycle_k_semicircularbeamdiscrete{found_it}.k; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('semicircularbeamdiscrete.m',params); %Load the equivalent discretized model. 
      k=assemble_system(net,'K',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      k=reordermatrixtop(k,node1ids,node2ids); %Reorder the matrix st.
      [k]=matrixcondensationtop(k); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_k_semicircularbeamdiscrete{found_it+1}.k=k; %Store k.
      recycle_k_semicircularbeamdiscrete{found_it+1}.w=params.w; %Store the w parameter for that k.
      recycle_k_semicircularbeamdiscrete{found_it+1}.h=params.h; %Store the h parameter for that k.
      recycle_k_semicircularbeamdiscrete{found_it+1}.radius=params.radius; %Store the radius parameter for that k.
      recycle_k_semicircularbeamdiscrete{found_it+1}.alpha=params.alpha; %Store the alpha parameter for that k.
   end
   %Rotate k into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*k*R')*To'; %Return the global 12*12 stiffness matrix.
case 'pos'
  % Compute relative positions of beam nodes
  x=params.radius*sin(params.alpha); %Relative x position of node2.
  y=params.radius*(1-cos(params.alpha)); %Relative y position of node2.
  [RL1,RL2]=rigidlinkcoordinates(params); %[x;y;z] coordinated of both rigidlinks wrt end nodes.
  P=RL2-RL1; %difference of rigidlink coordinates.
  output=R*[0,x+P(1);0,y+P(2);0,0+P(3)]; %R*[node1postions,node2postions]
case 'display'
   L=params.radius*sin(params.alpha); %projected curve length in local coordinates.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations. To.
   [RL1,RL2,r1,r2]=rigidlinkcoordinates(params); %[x;y;z] coordinated of both rigidlinks.
   q=To'*q;
   Rp=(nodes(1).pos-R*RL1);     
   displaycircbeam4(q, Rp, R, params.w, params.h, params.radius, params.alpha);   
case 'check'
   S=[];
   if ~isfield(params,'radius'),S=strcat(S,' Semicircularbeam model needs (radius).');end
   if ~isfield(params,'h'),S=strcat(S,' Semicircularbeam model needslayer thickness (h).');end
   if ~isfield(params,'w'),S=strcat(S,' Semicircularbeam model needs width (w).');end
   if ~isfield(params,'alpha'),S=strcat(S,' Semicircularbeam model needsangle sweep (alpha).');end
   if ~isfield(params,'Youngsmodulus'),S=strcat(S,' Semicircularbeam model needsprocess file Youngs modulus (Youngsmodulus).');end
   if ~isfield(params,'fluid'),S=strcat(S,' Semicircularbeam model needsprocess file viscous fluid layer thickness (fluid).');end
   if ~isfield(params,'viscosity'),S=strcat(S,' Semicircularbeam model needsprocess file viscosity of fluid layer (viscosity).');end
   if ~isfield(params,'density'),S=strcat(S,' Semicircularbeam model needsprocess file material density (density).');end
   output=[S];
otherwise
  output = [];
end


