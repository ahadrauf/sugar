%function macromatrix('mftriangle.m',{'A','B','C'},'K')
function [J]=macromatrix(netlist,netlistparameters,nodenames,matrixtype)
%by jvclark - July2001
net=cho_load(netlist,netlistparameters);
[M]=assemble_system(net,matrixtype);
dof=net.dof;
numnodes=length(nodenames);
numcoords=length(net.elements(2).node_info.dynamic{1,2})*3;
for i=1:numnodes %for each node
   G=zeros(dof+numcoords,dof+numcoords); %init MNA matrix
   G(1:dof,1:dof)=M; %enter matrix
   G((i-1)*numcoords+1:i*numcoords,dof+1:dof+numcoords)=eye(numcoords); %enter eye
   G(dof+1:dof+numcoords,(i-1)*numcoords+1:i*numcoords)=eye(numcoords); %enter eye
   for q=1:numcoords %for each coordinate
      Xs=zeros(dof+numcoords,1); %init source vector
      Xs(dof+q)=1; %set source displacement
      Qs(1:numcoords,1)=Xs(dof+1:dof+numcoords,1); %pick off displacement source vector
      V=G\Xs; %MNA solution vector
      F(1:numcoords,1)=V(dof+1:dof+numcoords,1); %pick off excitation vector
      for f=1:numcoords
         J(f,q)=-F(f)/Qs(q);
      end %f
   end %q
end %i


