function harm_sines2

w0=1; 
m=1; 
k=m*w0^2; 
c_critical=sqrt(4*k*m);
c=c_critical * 0.2;
gamma=c/2/m;
I=sqrt(-1);
wr=(w0^2-2*gamma^2)^0.5
F0=1;

w0,wr

w1=0.25;
w2=0.617;
w3=wr-(w0-wr);
w4=wr
w5=w0
w6=1.208;
w7=1.5;
w8=2.5;

w=w0

phi=atan(c*w/(k-m*w^2)); %phi=atan(2*gamma*w/(w0^2-w^2));
A=F0/m/((w0^2-w^2)^2 + 4*gamma^2*w^2)^0.5;

i=0;
for t=0:pi*3/w/1000:pi*3/wr
    i=i+1;
    X(i)=           (A*exp(-I*phi))*exp(I*w*t);%X(i)=(A)*exp(I*(w*t-phi));
    Fs(i)= -k*      (A*exp(-I*phi))*exp(I*w*t);
    Fd(i)= -c* I*w* (A*exp(-I*phi))*exp(I*w*t);
    F(i) = F0*                      exp(I*w*t);
    T(i)=t;
end
mX=max(real(X))
%figure(1);clf;
%plot(T,Fs,'b',T,Fd,'g',T,F,'r',T,X,'c',T,Fs+Fd+F,'k');grid on;

figure(1);clf;
%plot(T,Fs,'b',T,Fd,'g',T,F,'r',T,X,'c',T,Fd+F,'k');grid on;
plot(T,F,'r',T,X,'c');grid on;
Max=max(real(Fd+F))

%hold on;plot(4,5.9);plot(4,-5.9)

%T=0;X=0;i=0;
%for t=0:pi*3/w/10:pi*3/w
%    i=i+1;
%    X(i)=           (A*exp(-I*phi))*exp(I*w*t);
%    T(i)=t;
%end
%hold on;
%plot(T,X,'oc');
%phi
%grid on;

i=0;W=0;
for w=0:pi*3/wr/5000:pi*1/wr
    i=i+1;
    A(i)=F0/m/((w0^2-w^2)^2 + 4*gamma^2*w^2)^0.5;
    W(i)=w;
    v(i) = w*F0/(m*sqrt((w0^2-w^2)^2+4*gamma^2*w^2));
    dv(i)=F0/(m*sqrt((w0^2-w^2)^2+4*gamma^2*w^2))-1/2*w*F0*(-4*(w0^2-w^2)*w+8*gamma^2*w)/(m*((w0^2-w^2)^2+4*gamma^2*w^2)^(3/2));
end 
figure(2);plot(W,A,'b');grid on;
figure(6);plot(W,v,'b');grid on;
figure(7);plot(W,dv,'b');grid on;

%max(real(F+Fd))
%min(real(F+Fd))
%figure(3);plot(T,F+Fd);grid on;

i=0;W=0;
for w=0:0.01:5
i=i+1;
phase(i)=atan(2.*gamma.*w./(w0.^2-w.^2));
if phase(i)<0
    phase(i)=phase(i)+pi;
end
W(i)=w;
end
figure(3);
plot(W,phase);grid on;


