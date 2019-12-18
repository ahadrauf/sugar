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
