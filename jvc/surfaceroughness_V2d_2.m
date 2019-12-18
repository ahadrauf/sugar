%x_=0; y_=1; b_=10e-9; c_=10; A_=2e-6; V0_=15; V=surfaceroughness_V2d_2(x_,y_,b_,c_,A_,V0_)
%clear all; res=30; for i=0:res; y_(i+1)= pi*i/res; x_=0; ; b_=1; c_=3; A_=20e-6; V0_=15; V(i+1)=surfaceroughness_V2d_2(x_,y_(i+1),b_,c_,A_,V0_); end;figure(2);plot(V,y_);xlabel('V');ylabel('y');grid on;
%clear all; tic; res=40; for i=1:res; y_(i)= pi * i/(res+1); x_=0; b_=1; c_=10; A_=2e-6; V0_=15; [V(i),X(i),Y(i)]=surfaceroughness_V2d_2(2e-6,y_(i),b_,c_,A_,V0_); end; figure(1);plot(X,Y); grid on; toc
%figure(2);plot(V0_+b_*sin(c_*y_),y_);grid on; 

function [V,x,y]=surfaceroughness_V2d_2(x_,y_,b_,c_,A_,V0_,y1,y2)
syms A V0 c b x y n;
S = symsum(2*sinh(n*(A+x))*V0*sin(n*y)*(1+(-1)^(n+1))/(n*3.14159265358979*sinh(2*n*A)), n, 1, 2000) ;

b=b_; c=c_; V0=V0_; A=A_; y=y_; x=x_; 

V = S + sinh(c*(A+x))*b*sin(c*y)/sinh(2*c*A);

op=optimset('tolx',eps);
x=fzero('find_zero_surfroughness1',A,op,V-V0,b,c,V0,A,y);

X=x;
Y=y;
V=eval(V);



