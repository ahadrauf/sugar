%width=1e-6;gap=2e-6;z=-10e-6 + sqrt(-1)*gap;tol=1e-20;t1=-1;t2=-1e-50;t=sc7b(width,gap,z,tol,t1,t2)
function t=sc7b(width,gap,z,tol,t1,t2)
%function t=sc7b(width,gap,z,tol)
t=t1;

%gap type 2
g=gap;
p=width;
k=(1+p/g);
a=-1+2*k^2+ 2*k*sqrt(k^2-1); %for t=-1..0
%a=-1+2*k^2+ -2*k*sqrt(k^2-1) 
%R=sqrt((t+1)/(t+a));
%zoft='g/pi*((a+1)/sqrt(a)*atanh(R)+(a-1)/sqrt(a)*R/(1-R^2)+log((R*sqrt(a)-1)/(R*sqrt(a)+1))) - z'; %without R plugged in.

%t1=-eval(num2str(a,'%1.13f')); %<- slightly less than -a!
%t2=-1;
%t=t1;
zoft='g/pi*((a+1)/sqrt(a)*atanh(sqrt((t+1)/(t+a)))+(a-1)/sqrt(a)*sqrt((t+1)/(t+a))/(1-((t+1)/(t+a)))+log((sqrt((t+1)/(t+a))*sqrt(a)-1)/(sqrt((t+1)/(t+a))*sqrt(a)+1))) - z'; %gap type 2

%gap type 1
%g=gap;
%zoft='g/pi* (1 + t + log(t)) - z'; 

%set range
%t1=-1;t2=-1e-50;

maxiterations=1000;

if imag(eval(zoft)) > 0
   t1pos=1;
else
   t1pos=0;
end
i=0;
F=abs(imag(eval(zoft)));

if t1pos
   while F>tol & i<maxiterations
      t=(t2+t1)/2;
      if imag(eval(zoft)) > 0
         t1=t;
      else
         t2=t;
      end
      F=abs(imag(eval(zoft)));
      i=i+1;
   end
else
   while F>tol & i<maxiterations
      t=(t2+t1)/2;
      if imag(eval(zoft)) < 0
         t1=t;
      else
         t2=t;
      end
      F=abs(imag(eval(zoft)));
      i=i+1;
   end
end
iterations=i

zoftnoz='g/pi*((a+1)/sqrt(a)*atanh(sqrt((t+1)/(t+a)))+(a-1)/sqrt(a)*sqrt((t+1)/(t+a))/(1-((t+1)/(t+a)))+log((sqrt((t+1)/(t+a))*sqrt(a)-1)/(sqrt((t+1)/(t+a))*sqrt(a)+1)))'; 
Z=eval(zoftnoz)
