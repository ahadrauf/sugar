%clear all;i=0;for kmax=1:2:50; i=i+1;A=2;b=0.2;x=1;y=1;v(i)=surfaceroughness_V2d_1(10,A,b,x,y,kmax); end; clf;figure(1); plot(1:length(v),(v(length(v))-v)/v(length(v))  );hold on;grid on;

function V=surfaceroughness_V2d_1(V0,A,b,x,y,kmax)


V=0;
for k=1:kmax;
Vxy=(1/2*(-1+(-1)^k)*(-exp(k*(A+x))+exp(-k*(A+x)))*(2*V0+b*pi*k)*sin(k*y)/(pi*k*(exp(2*k*A)-exp(-2*k*A))));
V=V+Vxy;
end

