function [net,M,D,K,Gmna] = cho2mna(net)
%assemble G and extract the zero and nonzero parts
G=assemble_system(net,'G',0);
zeroindex=0;
nonzeroindex=0;
mdk=zeros(net.dof,1);
g=zeros(net.dof,1);

%compress zero and nonzero indeces
for i=1:net.dof
   nonzerorow=find(G(i,:)~=0);
   if nonzerorow
      nonzeroindex=nonzeroindex+1;
      g(net.dof-nonzeroindex+1)=i;
   else
      zeroindex=zeroindex+1;
      mdk(zeroindex)=i;
   end
end

%renumber the net ids
net2=net;
for j=1:length(net.elements)
   for i=1:net.dof
      if ~isempty(net.elements(j).var_ids)
         index=find(net.elements(j).var_ids==i);
         if ~isempty(index)
            if mdk(i)
               net2.elements(j).var_ids(index)=mdk(i);
            end
            if g(i)
               net2.elements(j).var_ids(index)=g(i);
            end
         end
      end
   end
end
for i=1:net.dof
   for j=1:length(net.nodes)
      coordnames=fieldnames(net.nodes(j).vars);
      for k=1:length(coordnames)      
         if getfield(net.nodes(j).vars,coordnames{k})==i         
            if mdk(i)
              net2.nodes(j).vars=setfield(net2.nodes(j).vars,coordnames{k},mdk(i));
            end
            if g(i)
               net2.nodes(j).vars=setfield(net2.nodes(j).vars,coordnames{k},g(i));
            end
         end
      end
   end
end

%set dof for 1st and 2nd order systems
net.dof2ndorder=length(find(mdk~=0));
net.dof1storder=length(find(g~=0));

%calculate new matrices using the new node and variable indices
k=assemble_system(net2,'K',0);
m=assemble_system(net2,'M',0);
d=assemble_system(net2,'D',0);
g=assemble_system(net2,'G',0);
K=k(1:net.dof2ndorder,1:net.dof2ndorder);
M=m(1:net.dof2ndorder,1:net.dof2ndorder);
D=d(1:net.dof2ndorder,1:net.dof2ndorder);
G=g(net.dof2ndorder+1:net.dof2ndorder+net.dof1storder,net.dof2ndorder+1:net.dof2ndorder+net.dof1storder);
Gmna = ...
   [M                                       D                                       zeros(net.dof2ndorder,net.dof1storder); 
    zeros(net.dof2ndorder,net.dof2ndorder)  eye(net.dof2ndorder)                    zeros(net.dof2ndorder,net.dof1storder);
    zeros(net.dof1storder,net.dof2ndorder)  zeros(net.dof1storder,net.dof2ndorder)  G];
 