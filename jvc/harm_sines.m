function harm_sines
figure(1);clf;

w0=1; m=1; 
k=m*w0^2; c=sqrt(4*k*m) *0.2;
%c=(1/2)*m*w0;

gamma=c/2/m;

i=0;
I=sqrt(-1);

w0
wr=(w0^2-2*gamma^2)^0.5
w=wr-(w0-wr);
%w=0.617;
%w=0.25;
%w=1.208;
%w=1.5;
%w=2.5;
%w=wr
%w=w0

x0=1;
F0=1;
phi=atan(c*w/(k-m*w^2));
%phi=atan(2*gamma*w/(w0^2-w^2));
A=F0/m/((w0^2-w^2)^2 + 4*gamma^2*w^2)^0.5;
for t=0:pi*3/w/1000:pi*3/w
    i=i+1;
    X(i)=           (A*exp(-I*phi))*exp(I*w*t);
    X(i)=           (A)*exp(I*(w*t-phi));
    Fs(i)= -k*      (A*exp(-I*phi))*exp(I*w*t);
    Fd(i)= -c* I*w* (A*exp(-I*phi))*exp(I*w*t);
    F(i) = F0*                      exp(I*w*t);
    T(i)=t;
end
mX=max(real(X))
plot(T,Fs,'b',T,Fd,'g',T,F,'r',T,X,'c',T,Fs+Fd+F,'k');
T=0;X=0;i=0;
for t=0:pi*3/w/10:pi*3/w
    i=i+1;
    X(i)=           (A*exp(-I*phi))*exp(I*w*t);
    T(i)=t;
end
hold on;
plot(T,X,'oc');
phi
grid on;

i=0;W=0;
for w=0:pi*3/w/5000:pi*4/w
    i=i+1;
    A(i)=F0/m/((w0^2-w^2)^2 + 4*gamma^2*w^2)^0.5;
    W(i)=w;
end 
%figure(2);
%plot(W,A,'b');
grid on;

