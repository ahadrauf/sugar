function harm_amp1
m=1000;
k=100;
c=100;
f0=10;
I=sqrt(-1);
w0=sqrt(k/m);
g=c/2/m;
wd=sqrt(w0^2-g^2);
wr=sqrt(w0^2-2*g*g);
w=wr;
phi=atan(2.*g.*w./(w0.^2-w.^2))
t=0:0.1:50;
A=f0./m./sqrt((w0.^2-w.^2).^2 + 4.*g.^2.*w.^2);
x = A.*exp(I.*(w.*t - phi));
v = A.*I.*w.*exp(sqrt(-1).*(w.*t - phi));
a = A.*((I.*w).^2).*exp(sqrt(-1).*(w.*t - phi));
figure(1);plot(t,x,'r',t,v,'b',t,a,'g');grid on;
figure(2);plot(t,imag(x),'r',t,imag(v),'b',t,imag(a),'g');grid on;

t=phi/w
x = A.*exp(I.*(w.*t - phi))
v = A.*I.*w.*exp(I.*(w.*t - phi))
a = A.*((I.*w).^2).*exp(I.*(w.*t - phi))
ma=m*a
cv=c*v
kx=k*x
F=f0*exp(I.*w.*t)
K=(F-ma-cv)/x



