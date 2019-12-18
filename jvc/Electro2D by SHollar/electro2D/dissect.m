function [array_of_segments] = dissectmex(line,delta)

%    function [array_of_segments] = dissect(line,delta)
%
%    dissect takes in a line of the form specified in electro2d and 
%    dissects it into an array of segments of length delta.  
%    This function is used by electro2d.
%
%    See also electro2d, plot_electro2d

distance = sqrt((line(2)-line(4))^2+(line(1)-line(3))^2);
num_seg = ceil(distance/delta);
angle = atan2(line(4)-line(2),line(3)-line(1));

for i=1:num_seg
   x_mean=line(1)+(line(3)-line(1))*(i-0.5)*delta/distance;
   y_mean=line(2)+(line(4)-line(2))*(i-0.5)*delta/distance;
   array_of_segments(i,:)=[x_mean y_mean angle line(5)];
end
