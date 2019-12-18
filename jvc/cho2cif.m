%Netlist to CIF converter
%Written by Andy Kuo

%the netlist to CIF converter is able to take a netlist after it is loaded into sugar.
%Added features include, etchholes, anchor, extend beam and display
%to call the function use the call cho2cif(net,{'p1'},0)  the net is the netlist that is loaded, p1 is the layer tha you  want
%and 0 is if you want extend beams on or off, in this case it is off.
function cho2cif(net,layer,ciffill) 
format compact 
figure(2);clf;hold on;grid on;axis equal; 
title('CIF'); 
xlabel('meters'); 
ylabel('meters'); 
%CIF headers
fid=fopen('cif.txt', 'w') 
fprintf(fid,'DS 1 1 1; \n')
fprintf(fid, '9 sugarcif; \n')
for lyr = 1 : length(layer)%runs through for layers
   fprintf(fid,'L POLY%d; \n',lyr)
   color= 1:length(layer)
   
   for i=1:length( net.elements ) %check all netlist elements 
      %begin polygon routine
      if isfield(net.elements(i).parameter, 'polygon')
         fprintf(fid,'P') %for the polygon
         numvertices=length(net.elements(i).parameter.polygon);
         for elm=1:numvertices
            ox=net.elements(i).parameter.ox;
            oy=net.elements(i).parameter.oy;
            oz=net.elements(i).parameter.oz;
            x=net.elements(i).parameter.polygon(elm,1);
            y=net.elements(i).parameter.polygon(elm,2);
            q=[x;y;0] + net.nodes(i).pos;
            p=rot2local(ox,oy,oz)'*q;
            fprintf(fid,' %d,%d', x, y)
         end
         fprintf(fid,';\n')
      end
      if isfield(net.elements(i).parameter,layer{lyr}) %must be right layer 
         L=net.elements(i).parameter.l*1e8; %get length 
         W=net.elements(i).parameter.w*1e8; %get width 
         theta=net.elements(i).parameter.oz; %get xy-plane orientation 
         id=net.elements(i).node_ids(1); %node index 
         X=net.nodes(id).pos(1)*1e8+L/2*cos(theta); %global X position 
         Y=net.nodes(id).pos(2)*1e8+L/2*sin(theta); %global Y position 
         a=cos(theta); 
         b=sin(theta);
         maxseparation=3000
         
         %Fill in beam corners 
         flagnode = {'C'}; %store nodes that are modified 
         for nodesource = 1:length( net.elements(i).node ) %get a node to check 
            for j = 1 : length(net.elements) %check netlist for another 
               if isfield(net.elements(j).parameter,layer{lyr}) %must be right layer 
                  for nodetest=1:length( net.elements(j).node) %get a 2nd node 
                     if (strcmp(net.elements(i).node{nodesource},net.elements(j).node{nodetest})&(i~=j)& isempty(strmatch(net.elements(i).node{nodesource}, flagnode,'exact'))) 
                        %has to be same node, but differnet element, and not already checked 
                        
                        flagnode{length(flagnode)+1}=net.elements(i).node{nodesource}; %flag it 
                        
                        if isfield(net.elements(j).parameter,'ciffill') 
                           netciffill = net.elements(j).parameter.ciffill; 
                           %override 
                           
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
                        end%netciffill   
                     end %strcmp 
                  end %nodetest 
               end %isfield
            end %j 
         end %nodesource 
         
         
         fprintf(fid,'B %d %d %d, %d %d,%d;\n',round(L),round(W),round(X),round(Y),round(cos(theta)*100),round(sin(theta)*100)); 
         %cif output
         B=rot2local(0,0,atan(b/a))*2\[[-L L L -L];[-W -W W W];[0 0 0 0]]; 
         patch((X+B(1,:))*1e-8,(Y+B(2,:))*1e-8, color(lyr)); 
         %display routine
         %if beam is larger than 40um on both sides, and larger than 60um on at least 1 side, etchholes will form
         if length(net.elements(i).node)==2 &(L>=6000 | W>=6000)%etchole
            if (L>4000&W>4000)
               fprintf(fid,'L HOLE1;\n')%places on right layer
               etchhole(maxseparation,L,W,fid,theta,X,Y)%call to the function etchhole
               fprintf(fid,'L POLY%d; \n',lyr)    %reverts back to old layer
            end
         end%etchole
         if length(net.elements(i).node)==1 %anchor
            
            margin=400%sets margin between poly and anchor to 4um
            l=abs(L-2*margin)
            w=abs(W-2*margin)
            fprintf(fid, 'L ANCHOR1; \n')%uses anchor layer
            D=rot2local(0,0,atan(b/a))*2\[[-l l l -l];[-w -w w w];[0 0 0 0]];
            patch((X+D(1,:))*1e-8,(Y+D(2,:))*1e-8, 4); 
            fprintf(fid,'B %d %d %d, %d %d,%d;\n',round(l),round(w),round(X),round(Y),round(cos(theta)*100),round(sin(theta)*100)); 
            
            fprintf(fid,'L POLY%d; \n',lyr)    
            
         end %anchor
         
         
      end %layer 
   end %i 
end
fprintf(fid,'DF; \n')
fclose(fid)