%harmonicosc_amp

m=1000;
k=100;
c=100;
F=10;
g=c/2/m;
w0=sqrt(k/m)
g=w0/sqrt(2) * 0.5;
%g=w0/sqrt(2) ;
%g=w0/sqrt(2) * 1.4;
wd=sqrt(w0^2-g^2)

wr=sqrt(w0^2-2*g*g)
w=0:3*w0/10000:3*w0;
A=F./m./sqrt((w0.^2-w.^2).^2 + 4.*g.^2.*w.^2);


figure(1);
plot(w,A,'b');
grid on;
hold on;


dA=-1/2.*F.*(-4.*w0.^2.*w+4.*w.^3+8.*g.^2.*w)./(m.*(w0.^4-2.*w0.^2.*w.^2+w.^4+4.*g.^2.*w.^2).^(3/2));
figure(2);
plot(w,dA,'b');
grid on;
%hold on;



x = A*cos(-w*t + phi);
v = A*sin(-w*t + phi)*w;
a = -A*cos(-w*t + phi)*w;

