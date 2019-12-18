function r=newton(z,t_guess,fnc,small)
width=2e-6;
gap=2e-6;
r=t_guess;
dx=1;
i=0;
while abs(dx) > small
   i=i+1;   
   if i>2000,r, error(' ::: EXCESSIVE ITERATIONS ::: '); end
   [dx,df,a]=feval(fnc,width,gap,r);
   dx=dx-z;
   if abs(df) < small, error(' ::: DIVERGING ::: '); end
   r=r+dx;
end
iterations=i