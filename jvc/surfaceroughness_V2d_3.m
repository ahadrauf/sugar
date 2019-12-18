%tic; b_=10e-9;c_=100;A_=2e-6;V0_=15; [V,x,y]=surfaceroughness_V2d_3(b_,c_,A_,V0_);toc, figure(1);plot(eval(x(2:10)),y(2:10),'-o'); grid on;

function [V,x,y]=surfaceroughness_V2d_3(b_,c_,A_,V0_)

%Solve for V(x,y) symbolically
for  number_summation_terms = [3 30 300 3000 30000]
tic;

    
    syms A V0 c b x y n;
Vxy_sym = symsum(2*sinh(n*(A+x))*V0*sin(n*y)*(1+(-1)^(n+1))/(n*3.14159265358979*sinh(2*n*A)), n, 1, number_summation_terms) + sinh(c*(A+x))*b*sin(c*y)/sinh(2*c*A);

b=b_; %perturbation in V0
c=c_; %number of cycles / 2
V0=V0_; %potential contour
A=A_; %gap / 2
%y=y_; %y coord
%x=x_; %x coord

if 0
i=0;
for y=0:pi/300:pi
    i=i+1;
    x=A;
    v11(i) = eval(Vxy_sym);
    y11(i)=y;
end
figure(1); plot(y11,v11);grid on;
end
if 1
%one period along length 
y=pi/2; %middle of the plate
L=pi; %length of the plate
lambda = 2*L/c; %length of one period
y1=y-lambda/2; %y min
y2=y+lambda/2; %y max

y1;
y2;
lambda;
c;
%Vxy_sym

res = 2;
i = 0;
op=optimset('tolx',eps);
for Y=y1:(y2-y1)/res:y2
    i = i + 1;

    %solve for x given y, V contour, and other params
    x22(i) = fzero('find_zero_surfroughness1', A, op, Vxy_sym-V0, b, c, V0, A, Y);
    y22(i) = Y;
toc    

end
%figure(2);plot(y22,x22,'-ob');grid on;

xmin=abs(min(x22) - A)
xmax=abs(max(x22) - A)

%X=x;
%Y=y;
%V=eval(V);

end
end