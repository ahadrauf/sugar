function expnode2(px,py,height)
fprintf('\nWait while I plot data...\n\n');
%make a plot of the *.mat experimental data
figure(1);%surface(px,py,real(height));
fprintf('Then click on the figure at 2 distinct points\n');
fprintf('and I will plot the best fitting line bewteen those points.\n\n');
X=px;
Y=py;
Z=real(height);
%mouse input
[gx,gy]=ginput(2);
[leny,lenx]=size(Z);
%set up index bounds
i1=round(lenx/abs(max(max(X))-min(min(X)))*gx(1)+1-lenx/abs(max(max(X))-min(min(X)))*min(min(X)));
j1=round(leny/abs(max(max(Y))-min(min(Y)))*gy(1)+1-leny/abs(max(max(Y))-min(min(Y)))*min(min(Y)));
i2=round(lenx/abs(max(max(X))-min(min(X)))*gx(2)+1-lenx/abs(max(max(X))-min(min(X)))*min(min(X)));
j2=round(leny/abs(max(max(Y))-min(min(Y)))*gy(2)+1-leny/abs(max(max(Y))-min(min(Y)))*min(min(Y)));
resolution=max(abs(i2-i1),abs(j2-j1));
%increments along straight line
inci=sign(i2-i1)*abs(i2-i1)/resolution;
incj=sign(j2-j1)*abs(j2-j1)/resolution;
%get polynomial data
ii=0;
for k=1:resolution+1
   ii=ii+1;
   i=round(i1+(k-1)*inci);
   j=round(j1+(k-1)*incj);
   z(ii)=Z(j,i);   
   if ii>
   x(ii)=sqrt((X(j,i)-X(j1,i1))^2+(Y(j,i)-Y(j1,i1))^2);
end   
%fit data
p=polyfit(x,z,1);
%plot data contour
figure(2);clf;
plot(x,z,'b');
%plot fit contour
hold on;
grid on;
j=p(1)*x(1:resolution+1)+p(2);
plot(x(1:resolution+1),j,'r');
hold off;
%printout results
fprintf('      actual 1st mouse coords is (%f,%f)\n',gx(1),gy(1));
fprintf('interpolated 1st mouse coords is (%f,%f)\n',X(j1,i1),Y(j1,i1));
fprintf('      actual 2nd mouse coords is (%f,%f)\n',gx(2),gy(2));
fprintf('interpolated 2nd mouse coords is (%f,%f)\n',X(j2,i2),Y(j2,i2));
fprintf('       Best fit straight line is %f * X = %f\n',p(1),-p(2));
fprintf('                        Slope is %f\n',p(1));
fprintf('                        Angle is %f rad\n',atan(p(1)));

