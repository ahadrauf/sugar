function xx7
%x_shapefunctionscircbeam

R=100;
alpha=pi/4
L=R*sin(alpha);

q1=0; %y1
q2=0; %oz1
q2=tan(q2); %oz1
q3=0; %y2
q4=0; %oz2
q4=tan(q4); %oz2

i=0;
for x=0:L/100:L
   i=i+1;
   psi1=1-3*(x/L)^2+2*(x/L)^3;
   psi2=x*(1-x/L)^2;
   psi3=3*(x/L)^2-2*(x/L)^3;
   psi4=x^2/L*(x/L-1);
   X(i)=x;
   Y1(i) = psi1*q1 + psi2*q2 + psi3*q3 + psi4*q4;
   Y2(i)=(R-sqrt(R^2-x^2))+Y1(i);
%   Y2(i)=(R-sqrt(R^2-x^2));
end

figure(1);
plot(X,Y1,'b',X,Y2,'r');
%plot(X,Y1,'b');
grid on;
axis equal;

s=linspace(0,L,20);
[cb]=circularshape(R,alpha,0,0,0,0,s)
hold on
plot(s,cb,'g');

function [cb]=circularshape(R,alpha,y1,oz1,y2,oz2,x)
L=R*sin(alpha);
q1=y1;
q2=oz1;
q2=tan(q2); %angle to slope
q3=y2;
q4=oz2;
q4=tan(q4); %angle to slope
psi1=1-3*(x./L).^2+2*(x./L).^3;
psi2=x.*(1-x./L).^2;
psi3=3*(x./L).^2-2*(x./L).^3;
psi4=x.^2./L.*(x./L-1);
Y = psi1.*q1 + psi2.*q2 + psi3.*q3 + psi4.*q4;
cb=(R-sqrt(R.^2-x.^2))+Y;
