function cho2cif3(net,layer,ciffill)
format compact
figure(2);clf;hold on;grid on;axis equal;
title('CIF');
xlabel('meters');
ylabel('meters');
fid=fopen('cif.txt', 'w')
for i=1:length( net.elements ) %check all netlist elements
   if isfield(net.elements(i).parameter,layer) %must be right layer
      L=net.elements(i).parameter.l*1e8; %get length
      W=net.elements(i).parameter.w*1e8; %get width
      theta=net.elements(i).parameter.oz; %get xy-plane orientation
      id=net.elements(i).node_ids(1); %node index
      X=net.nodes(id).pos(1)*1e8+L/2*cos(theta); %global X position
      Y=net.nodes(id).pos(2)*1e8+L/2*sin(theta); %global Y position
      a=cos(theta);
      b=sin(theta);
      %Fill in beam corners
      flagnode = {'C'}; %store nodes that are modified
      for nodesource = 1:length( net.elements(i).node ) %get a node to check
         for j = 1 : length(net.elements) %check netlist for another
            if isfield(net.elements(j).parameter,layer) %must be right layer
               for nodetest=1:length( net.elements(j).node) %get a 2nd node
                  if (strcmp(net.elements(i).node{nodesource},net.elements(j).node{nodetest})&(i~=j) & isempty(strmatch(net.elements(i).node{nodesource}, flagnode,'exact'))) %has to be same node, but differnet element, and not already checked
                     flagnode{length(flagnode)+1}=net.elements(i).node{nodesource}; %flag it
                     
                     if isfield(net.elements(j).parameter,'ciffill')
                        netciffill = net.elements(j).parameter.ciffill; %override
                        
                     else
                        
                        netciffill = ciffill; %no override
                        
                     end
                     
                     if netciffill
                        W2=net.elements(j).parameter.w *1e+08; %get width of 2nd element
                        
                        if (nodesource == 1)
                           
                           
                           L=L+W2/2; %adjust lenght of 1st element
                           X=X-((W2/4)*cos(theta)); %adjust X position of 1st element
                           Y=Y-((W2/4)*sin(theta)); %adjust Y position of 1st element
                        end
                        if (nodesource==2)
                           L=L+W2/2; %adjust lenght of 1st element
                           X=X+((W2/4)*cos(theta)); %adjust X position of 1st element
                           Y=Y+((W2/4)*sin(theta)); %adjust Y position of 1st element
                        end
                        
                        
                     end %strcmp
                  end %nodetest
               end %layer
            end %j
         end %nodesource
      end %corners
      
      fprintf(fid,'B %d %d %d, %d %d,%d;\n',round(L),round(W),round(X),round(Y),round(cos(theta)*100),round(sin(theta)*100)); %cif
      B=rot2local(0,0,atan(b/a))*2\[[-L L L -L];[-W -W W W];[0 0 0 0]];
      patch((X+B(1,:))*1e-8,(Y+B(2,:))*1e-8,3);
      
      
   end %layer
end %i
fclose(fid)
