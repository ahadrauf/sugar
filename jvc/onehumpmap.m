function [LAMBDA,XINIT,FLIMIT]=onehumpmap

i = 0;
for lambda = 0 : 0.01 : 4
   i = i + 1;
   j = 0;
   for initial_x = -1 : 0.01 : 2
      x = initial_x;
      j = j + 1;
      for iteration = 1 : 100
         x = lambda*x*(1-x);         
      end
      if x < -1
         x = -1;
      end
      FLIMIT(i,j) = x; 
      LAMBDA(i,j) = lambda;
      XINIT(i,j) = initial_x;
   end
end

figure(1);
clf;
%hold on;
%surf(LAMBDA,XINIT,FLIMIT);
%shading interp;
colormap default;
%mesh(LAMBDA,XINIT,FLIMIT);
contourf(LAMBDA,XINIT,FLIMIT);
xlabel('lambda');
ylabel('x init');
zlabel('F limit');
%hold off;
