function shownodes(net)

for i=1:length(net.nodes)
   x=net.nodes(i).pos;
   name=net.nodes(i).name;
   text(x(1),x(2),x(3),name);
end

   