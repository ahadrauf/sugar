
function etchholecoords(maxseparation,L,W,theta,X,Y)
%maxseparation is maximum distance between etch holes
%L is the length of the big box
%W is the width of the big box
%theta is the oz angle
%X is the middle x-coordinate of the big box
%Y is the middle y-coordinate of the big box
color=3; %display color
minetchholesize=5e-6; %minimum allowable etch hole dimension
l=minetchholesize; %etch hole length
w=minetchholesize; %etch hole width
B=rot2local(0,0,theta)*2\[[-L L L -L];[-W -W W W];[0 0 0 0]];
patch(X+B(1,:),Y+B(2,:),color); %show big box
for xindex=1:ceil(L/maxseparation)-1
   x(xindex)=(L/ceil(L/maxseparation)*xindex); 
   for yindex=1:ceil(W/maxseparation)-1
      y(yindex)=((W/ceil(W/maxseparation)*yindex)-W/2); 
      B=rot2local(0,0,theta)*2\[[-l l l -l];[-w -w w w];[0 0 0 0]];
      p=rot2local(0,0,theta)'*[x(xindex)-L/2;y(yindex);0]; %point
		patch(p(1)+X+B(1,:),p(2)+Y+B(2,:),color); %display
   end
end
