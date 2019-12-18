function r=secant(x,dx)

r=x;
dx0=dx;
dx1=dx0;

f0=feval('fnc',r);
i=0;
while abs(dx1)>1e-6
   r=r+dx0;
   f1=feval('fnc',r);
   dx1=-dx0*f1/(f1-f0);
   dx0=dx1;
   f0=f1;
   i=i+1;
end
iterations=i


function r=secant(x,dx)

r=x;
dx0=dx;
dx1=dx0;

f0=feval('fnc',r);
i=0;
while abs(dx1)>1e-6
   r=r+dx0;
   f1=feval('fnc',r);
   dx1=-dx0*f1/(f1-f0);
   dx0=dx1;
   f0=f1;
   i=i+1;
end
iterations=i

function f=fnc(r)
x=r;
S=3;alpha=0.7;Z=4i;
zoft='S*(x^(1-alpha)/(1-alpha)+x^(-alpha)/alpha)-Z';
%zoft=strrep(zoft,'S',num2str(S));
%zoft=strrep(zoft,'alpha',num2str(alpha));
%zoft=strrep(zoft,'Z',num2str(Z));
%zoft=strrep(zoft,'x',num2str(x));
%f=inline(zoft);t=fzero(f,5)
f=eval(zoft);
