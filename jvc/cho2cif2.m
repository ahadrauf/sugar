%net=cho_load('cant2.m');
%cho2cif(net,'p1');

function cho2cif(net,layer)
for i=1:length( net.elements )
   if isfield(net.elements(i).parameter,layer)
      L=net.elements(i).parameter.l*1e8;
      W=net.elements(i).parameter.w*1e8;
      O=net.elements(i).parameter.oz;
      id=net.elements(i).node_ids(1);
      X=net.nodes(id).pos(1)*1e8+L/2*cos(O);
      Y=net.nodes(id).pos(2)*1e8+L/2*sin(O);
      fprintf('B %d %d %d, %d %d, %d;\n',round(L),round(W),round(X),round(Y),round(cos(O)*100),round(sin(O)*100));
   end
end 