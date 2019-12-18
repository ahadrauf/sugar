

clear all functions;
clear all;
m=moviein(10);
tic
T =[...
 -3.2252e+002 -2.2904e+002  4.1011e-015  2.7578e+002;
  8.0046e+001 -7.0116e+001  1.4624e+002 -7.8085e+001;
  2.7915e+002 -2.4452e+002 -4.1934e+001  2.8701e+003;
            0            0            0  1.0000e+000];
view(T);
         
for i=1:10
   p.f=(i)*10e-6;
   clf;
   net=cho_load('robot3.m',p);
   q=cho_dc(net);
   cho_display2(net,q);
   view(-T);
   axis=([-6e-4 6e-4,  -5e-4 5e-4, -4e-4 4e-4 ]);
   axis equal;
   axis off;
   drawnow
   m(:,i)=getframe;
end
toc
%movie(m);
