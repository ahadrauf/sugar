function t=sc7(width,gap,z,tol)

%gap type 2
%g=gap;
%p=width;
%k=1+p/g;
%a=-1+2*k^2+ 2*k*sqrt(k^2-1);
%R=sqrt((t+1)/(t+a));
%zoft='g/pi*((a+1)/sqrt(a)*atanh(R)+(a-1)/sqrt(a)*R/(1-R^2)+log((R*sqrt(a)-1)/(R*sqrt(a)+1))) - z'; %without R plugged in.
%zoft='g/pi*((a+1)/sqrt(a)*atanh(sqrt((t+1)/(t+a)))+(a-1)/sqrt(a)*sqrt((t+1)/(t+a))/(1-((t+1)/(t+a)))+log((sqrt((t+1)/(t+a))*sqrt(a)-1)/(sqrt((t+1)/(t+a))*sqrt(a)+1))) - z'; %gap type 2

%gap type 1
g=gap;
zoft='g/pi* (1 + t + log(t)) - z'; 

%set range
t1=-1;t2=-1e-50;

t=t1;   
maxiterations=1000;

if real(eval(zoft)) > 0
   t1pos=1;
else
   t1pos=0;
end
i=0;
F=abs(real(eval(zoft)));
if t1pos
   while F>tol & i<maxiterations
      t=(t2+t1)/2;
      if real(eval(zoft)) > 0
         t1=t;
      else
         t2=t;
      end
      F=abs(real(eval(zoft)));
      i=i+1;
   end
else
   while F>tol & i<maxiterations
      t=(t2+t1)/2;
      if real(eval(zoft)) < 0
         t1=t;
      else
         t2=t;
      end
      F=abs(real(eval(zoft)));
      i=i+1;
   end
end
%iterations=i

