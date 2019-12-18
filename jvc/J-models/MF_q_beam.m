%net=cho_load('predisp1.m')
function [output] = MF_q_beam(flag, R, params, q, t, nodes, varargin)


switch(flag)
case 'vars'
   output.dynamic={ 1 {'x','y','z','rx','ry','rz'}; 2 {'x','y','z','rx','ry','rz'}}; 
case 'M' 
   global recycle_m_predisplacedbeamdiscrete
   num_m=length(recycle_m_predisplacedbeamdiscrete); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_m%Go through each one
      if (recycle_m_predisplacedbeamdiscrete{i}.w==params.w & ... %Check width.
         recycle_m_predisplacedbeamdiscrete{i}.qox1==params.qox1 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.qoy1==params.qoy1 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.qoz1==params.qoz1 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.qox2==params.qox2 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.qoy2==params.qoy2 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.qoz2==params.qoz2 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.qx2==params.qx2 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.qy2==params.qy2 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.qz2==params.qz2 & ... 
         recycle_m_predisplacedbeamdiscrete{i}.l==params.l & ... 
         recycle_m_predisplacedbeamdiscrete{i}.h==params.h) 
         found_it=i; %Store the cell index of the right one found.
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated k.
      m=recycle_m_predisplacedbeamdiscrete{found_it}.m; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('predisplacedbeamdiscrete.m',params); %Load the equivalent discretized model. 
      m=assemble_system(net,'M',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      m=reordermatrixtop(m,node1ids,node2ids); %Reorder the matrix st.
      [m]=matrixcondensationtop(m); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qox1=params.qox1;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qoy1=params.qoy1;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qoz1=params.qoz1;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qox2=params.qox2;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qoy2=params.qoy2;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qoz2=params.qoz2;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qx2=params.qx2;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qy2=params.qy2;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.qz2=params.qz2;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.l=params.l;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.m=m;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.w=params.w;
      recycle_m_predisplacedbeamdiscrete{found_it+1}.h=params.h;
   end
   %Rotate k into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*m*R')*To'; %Return the global 12*12 stiffness matrix.
case 'D'
   global recycle_d_predisplacedbeamdiscrete
   num_d=length(recycle_d_predisplacedbeamdiscrete); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_d%Go through each one
      if (recycle_d_predisplacedbeamdiscrete{i}.w==params.w & ... %Check width.
         recycle_d_predisplacedbeamdiscrete{i}.qox1==params.qox1 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.qoy1==params.qoy1 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.qoz1==params.qoz1 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.qox2==params.qox2 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.qoy2==params.qoy2 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.qoz2==params.qoz2 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.qx2==params.qx2 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.qy2==params.qy2 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.qz2==params.qz2 & ... 
         recycle_d_predisplacedbeamdiscrete{i}.l==params.l & ... 
         recycle_d_predisplacedbeamdiscrete{i}.h==params.h) 
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated d.
      d=recycle_d_predisplacedbeamdiscrete{found_it}.d; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('predisplacedbeamdiscrete.m',params); %Load the equivalent discretized model. 
      d=assemble_system(net,'D',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      d=reordermatrixtop(d,node1ids,node2ids); %Reorder the matrix st.
      [d]=matrixcondensationtop(d); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qox1=params.qox1;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qoy1=params.qoy1;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qoz1=params.qoz1;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qox2=params.qox2;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qoy2=params.qoy2;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qoz2=params.qoz2;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qx2=params.qx2;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qy2=params.qy2;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.qz2=params.qz2;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.l=params.l;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.h=params.h;      
      recycle_d_predisplacedbeamdiscrete{found_it+1}.d=d;
      recycle_d_predisplacedbeamdiscrete{found_it+1}.w=params.w;
   end
   %Rotate k into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*d*R')*To'; %Return the global 12*12 stiffness matrix.
case 'K'

net = varargin{1};
elt_no = varargin{2};
P = q_beam_parameters(net,elt_no);

%net.nodes(3)
% name: 'contact2'
% elt_ids: [3 4]
% pos: [3x1 double]
    
%net.elements(3)
% name: 'contact.anon1'
% model: 'MF_q_beam'
% node_ids: [2 3]
% parameter: [1x1 struct]
% R: [3x3 double]
% var_ids: [1 2 3 4 5 6 7 8 9 10 11 12]




    
    
    
    
   global recycle_k_predisplacedbeamdiscrete
   num_k=length(recycle_k_predisplacedbeamdiscrete); %Find the number of these used
   found_it=0; %Initialize. The right stiffness is not found yet.
   for i=1:num_k%Go through each one
      if (recycle_k_predisplacedbeamdiscrete{i}.w==params.w & ... %Check width.
         recycle_k_predisplacedbeamdiscrete{i}.qox1==params.qox1 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.qoy1==params.qoy1 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.qoz1==params.qoz1 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.qox2==params.qox2 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.qoy2==params.qoy2 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.qoz2==params.qoz2 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.qx2==params.qx2 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.qy2==params.qy2 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.qz2==params.qz2 & ... 
         recycle_k_predisplacedbeamdiscrete{i}.l==params.l & ... 
         recycle_k_predisplacedbeamdiscrete{i}.h==params.h) 
         found_it=i; %Store the cell index of the right one found.
         break; %Get out of the loop.
      end
   end
   if found_it %Then just use the precalculated k.
      k=recycle_k_predisplacedbeamdiscrete{found_it}.k; %This only checks geometry for now. Add process checking later!
   else %Then this k hasn't been made before, so make it from scratch and store it globally.
      net=cho_load('predisplacedbeamdiscrete.m',params); %Load the equivalent discretized model. 
      k=assemble_system(net,'K',1); %Assemble the large DOF matrix.
      node1ids=[lookup_coord(net,'node1','x'),lookup_coord(net,'node1','y'),lookup_coord(net,'node1','z'),lookup_coord(net,'node1','rx'),lookup_coord(net,'node1','ry'),lookup_coord(net,'node1','rz')]; %Find the matrix IDs for first node.
      node2ids=[lookup_coord(net,'node2','x'),lookup_coord(net,'node2','y'),lookup_coord(net,'node2','z'),lookup_coord(net,'node2','rx'),lookup_coord(net,'node2','ry'),lookup_coord(net,'node2','rz')]; %Find the matrix IDs for last node.
      k=reordermatrixtop(k,node1ids,node2ids); %Reorder the matrix st.
      [k]=matrixcondensationtop(k); %Transform k*q to khat*q1, an effective 12*12 matrix.
      %Store it globally.
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qox1=params.qox1;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qoy1=params.qoy1;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qoz1=params.qoz1;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qox2=params.qox2;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qoy2=params.qoy2;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qoz2=params.qoz2;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qx2=params.qx2;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qy2=params.qy2;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.qz2=params.qz2;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.k=k;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.l=params.l;
      recycle_k_predisplacedbeamdiscrete{found_it+1}.h=params.h;      
      recycle_k_predisplacedbeamdiscrete{found_it+1}.w=params.w;
   end
   %Rotate k into the global frame.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   R=[R,o3;o3,R]; R=[R,o6;o6,R]; %12*12 rotation to global.   
   output=To*(R*k*R')*To'; %Return the global 12*12 stiffness matrix.
case 'pos'
   params.qx1=0;
   params.qy1=0;
   params.qz1=0;
   L=params.l;
   %Hermite polynomial coefficients
   yf0=params.qy1;
   yf00=params.qoz1;
   yfL=params.qy2;
   yfLL=params.qoz2;
   yf0L=(yfL-yf0)/L;
   yf00L=(yf0L-yf00)/L;
   yf0LL=(yfLL-yf0L)/L;
   yf00LL=(yf0LL-yf00L)/L;
   zf0=params.qz1;
   zf00=params.qoy1;
   zfL=-params.qz2;
   zfLL=-params.qoy2;
   zf0L=(zfL-zf0)/L;
   zf00L=(zf0L-zf00)/L;
   zf0LL=(zfLL-zf0L)/L;
   zf00LL=(zf0LL-zf00L)/L;
   s=L;
   hx=params.qx1+(1+(params.qx2-params.qx1)/L)*s;
   hy=yf0+s*(yf00+s*(yf00L+(s-L)*yf00LL));
   hz=zf0+s*(zf00+s*(zf00L+(s-L)*zf00LL));
   % Compute relative positions of beam nodes
   [RL1,RL2]=rigidlinkcoordinates(params); %[x;y;z] coordinated of both rigidlinks.
   P=RL2-RL1; %difference of rigidlink coordinates.
   output=R*[0,hx+P(1);0,hy+P(2);0,hz+P(3)]; %R*[node1postions,node2postions]
case 'display'
   [RL1,RL2,r1,r2]=rigidlinkcoordinates(params); %[x;y;z] coordinated of both rigidlinks.
   To=rigidlinkmatrix(params,R); %Compute the rigid link transformations. To.
   q=To'*q; %Transition to end nodes.
   Rp=(nodes(1).pos-R*RL1); %Transition to node1 coord.
   o3=zeros(3); o6=zeros(6); %Make the zero placeholders.
   r=[R,o3;o3,R]; r=[r,o6;o6,r]; %12*12 rotation to global.   
   Q=r*[0;0;0;params.qox1;params.qoy1;params.qoz1;params.qx2;params.qy2;params.qz2;params.qox2;params.qoy2;params.qoz2]; %Predisplaced picture.
   q=q+Q;
   displaybeam(q,Rp,R,params.l,params.w,params.h);   
case 'check'
   S=[];
   if ~isfield(params,'h'),S=strcat(S,' predisplacedbeamdiscrete model needslayer thickness (h).');end
   if ~isfield(params,'w'),S=strcat(S,' predisplacedbeamdiscrete model needs width (w).');end
   if ~isfield(params,'l'),S=strcat(S,' predisplacedbeamdiscrete model needsangle sweep (alpha).');end
   if ~isfield(params,'Youngsmodulus'),S=strcat(S,' predisplacedbeamdiscrete model needsprocess file Youngs modulus (Youngsmodulus).');end
   if ~isfield(params,'fluid'),S=strcat(S,' predisplacedbeamdiscrete model needsprocess file viscous fluid layer thickness (fluid).');end
   if ~isfield(params,'viscosity'),S=strcat(S,' predisplacedbeamdiscrete model needsprocess file viscosity of fluid layer (viscosity).');end
   if ~isfield(params,'density'),S=strcat(S,' predisplacedbeamdiscrete model needsprocess file material density (density).');end
   output=[S];
otherwise
  output = [];
end

