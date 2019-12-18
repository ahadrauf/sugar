function [fmax]=onehumpmap3(x,lambda,nmax)

%lambda=3;
figure(2);
clf;
hold on;

f(1) = lambda*x*(1-x);         
for n = 2 : nmax
   f(n) = lambda*f(n-1)*(1-f(n-1));         
end
plot([1:nmax],f);
fmax=f(nmax);
avef=sum(f)/nmax
hold off;
