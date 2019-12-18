function plot_electro2d(seg1,seg2,delta_approx);

%    function plot_electro2d(seg1,seg2,delta_approx);
%
%    plot_electro2d plots the approximation of seg1 
%    and seg2.  Each line segment in seg1 and seg2 is
%    approximated as an array of line segments of 
%    length delta_approx. 
%   
%    see also electro2d


% clear original plot figure
clf;

% break up segments into smaller segments of length delta_approx
adding_matrix=[0 0 0 0];
for	n=1:size(seg1,1)
   adding_matrix=[adding_matrix;dissect(seg1(n,:),delta_approx)];
end
volt_m1=adding_matrix(2:size(adding_matrix,1),:);

adding_matrix=[0 0 0 0];
for	n=1:size(seg2,1)
   adding_matrix=[adding_matrix;dissect(seg2(n,:),delta_approx)];
end
volt_m2=adding_matrix(2:size(adding_matrix,1),:);

size_volt1 = size(volt_m1,1);
size_volt2 = size(volt_m2,1);

% total number of segments of length delta_approx
total_segment_number=size_volt1+size_volt2


% plot all lines
length=delta_approx/2;
plot_matrix=volt_m1;

for n=1:size(plot_matrix,1);
	theta=plot_matrix(n,3);
	x_c=plot_matrix(n,1);
	y_c= plot_matrix(n,2);
	line(length*cos(theta)*[1 -1] + x_c, ...
   	length*sin(theta)*[1 -1] + y_c);
end

plot_matrix=volt_m2;

for n=1:size(plot_matrix,1);
	theta=plot_matrix(n,3);
	x_c=plot_matrix(n,1);
	y_c= plot_matrix(n,2);
	line(length*cos(theta)*[1 -1] + x_c, ...
   	length*sin(theta)*[1 -1] + y_c);
end

   











