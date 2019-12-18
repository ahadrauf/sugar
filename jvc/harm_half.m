function harm_half
F0=1;
m=1;
gamma=0.1;
w0=1;
w=0:0.0001:3;
Amax=F0./m./(2*gamma.*sqrt(w0.^2-gamma.^2));
A=Amax.*gamma./sqrt((w0-w).^2.+gamma.^2);
Atrue=(F0./m)./sqrt((w0.^2-w.^2).^2.+4.*w.^2.*gamma.^2);

figure(1);clf;
plot(w,A,'b');grid on;hold on;
plot(w,Atrue,'g');

w1=w0+gamma;
A=Amax.*gamma./sqrt((w0-w1).^2.+gamma.^2);
plot(w1,A,'ro');
w2=w0-gamma;
A=Amax.*gamma./sqrt((w0-w2).^2.+gamma.^2);
plot(w2,A,'ro');

A
A^2