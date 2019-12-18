function harm_sines
m=1000;
k=100;
c=100;
F=10;
g=c/2/m;
w0=sqrt(k/m)
g=w0/sqrt(2) * 0.5;

figure(1);clf;
i=0;
x0=1;
I=sqrt(-1);
w=w0;
F0=1;
phi=atan(c*w/(k-m*w^2));
for t=0:0.01:7
    i=i+1;
    Fs(i)=-k*(x0*exp(-I*phi))*exp(I*w*t);
    Fd(i)=-c*I*w*(x0*exp(-I*phi))*exp(I*w*t);
    F(i)=F0*exp(I*w*t);
    X(i)=(x0*exp(-I*phi))*exp(I*w*t);
    T(i)=t;
end
plot(T,Fs,'b',T,Fd,'g',T,F,'r',T,X,'m');
phi
grid on;

figure(2);
plot(T,Fs+Fd+F,'m');
grid on;

mX=max(real(X))