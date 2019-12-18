function t=findroot1(fn,t1,t2,tol)
t=t1;   
maxiterations=1000;
if real(eval(fn)) > 0
   t1pos=1;
else
   t1pos=0;
end
i=0;
F=abs(real(eval(fn)));
if t1pos
   while F>tol & i<maxiterations
      t=(t2+t1)/2;
      if real(eval(fn)) > 0
         t1=t;
      else
         t2=t;
      end
      F=abs(real(eval(fn)));
      i=i+1;
   end
else
   while F>tol & i<maxiterations
      t=(t2+t1)/2;
      if real(eval(fn)) < 0
         t1=t;
      else
         t2=t;
      end
      F=abs(real(eval(fn)));
      i=i+1;
   end
end
iterations=i