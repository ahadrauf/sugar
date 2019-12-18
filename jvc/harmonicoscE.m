% harmonicoscE

%
if (1)
figure(3);clf;hold on;grid on;
m=1000;
k=100;
c=sqrt(4*m*k).*4;
m1=m;k1=k;c1=c;w0=sqrt(k/m);
t=0:0.25/4:30;
A2=0;A1=1;
gamma1=-(-c/2/m+(c*c/4/m/m - k/m)^0.5);
gamma2=-(-c/2/m-(c*c/4/m/m - k/m)^0.5);

xod=A1.*exp(-gamma1.*t)+A2.*exp(-gamma2.*t);
vod = -A1.*gamma1.*exp(-gamma1.*t)-A2.*gamma2.*exp(-gamma2.*t);
E = (1/2).*m.*vod.^2 + (1/2).*k.*xod.^2;
dEdt = -c.*vod.^2;
plot(t,xod,'b');
plot(t,vod./abs(vod(1)),'r');
plot(t,E./E(1),'g');
plot(t,dEdt./abs(dEdt(1)),'c');

m=m1*10;k=k1;c=c1;
gamma1=-(-c/2/m+(c*c/4/m/m - k/m)^0.5);
gamma2=-(-c/2/m-(c*c/4/m/m - k/m)^0.5);
xod=A1.*exp(-gamma1.*t)+A2.*exp(-gamma2.*t);
%plot(t,xod,'g');
m=m1;k=k1;c=c1*10;
gamma1=-(-c/2/m+(c*c/4/m/m - k/m)^0.5);
gamma2=-(-c/2/m-(c*c/4/m/m - k/m)^0.5);
xod=A1.*exp(-gamma1.*t)+A2.*exp(-gamma2.*t);
%plot(t,xod,'c');
m=m1;k=k1*10;c=c1;
gamma1=-(-c/2/m+(c*c/4/m/m - k/m)^0.5);
gamma2=-(-c/2/m-(c*c/4/m/m - k/m)^0.5);
xod=A1.*exp(-gamma1.*t)+A2.*exp(-gamma2.*t);
%plot(t,xod,'r');
end

if (1)
figure(2);clf;hold on;grid on;
t=0:0.25/4:30;
m=1000;
k=100;
c=100;
k1=k;m1=m;c=sqrt(4*m*k);A=0;B=1;gamma=c/2/m;

xcd=(A.*t + B).*exp(-gamma.*t);
vcd=A.*exp(-gamma.*t)-(A.*t+B).*gamma.*exp(-gamma1.*t);
E = (1/2).*m.*vcd.^2 + (1/2).*k.*xcd.^2;
dEdt = -c.*vcd.^2;
plot(t,xcd,'b');
plot(t,vcd./abs(vcd(1)),'r');
plot(t,E./abs(E(1)),'g');
plot(t,dEdt/abs(dEdt(1)),'c');

m=m1*10;k=k1;c=sqrt(4*m*k);
A=0;B=1;gamma=c/2/m;
xcd=(A.*t + B).*exp(-gamma.*t);
%plot(t,xcd,'g');grid on;
m=m1;k=k1*10;c=sqrt(4*m*k);
A=0;B=1;gamma=c/2/m;
%xcd=(A.*t + B).*exp(-gamma.*t);
%plot(t,xcd,'r');grid on;
end

if (0)
figure(1);clf;hold on;grid on; 
t=0:0.25/4:300;
m=10000;
k=100;
c=100;
m1=m;k1=k;c1=c;w0=sqrt(k/m);A=1;gamma=c/2/m;
wd=(w0^2-gamma^2)^0.5;
xud = exp(-gamma.*t) .* A .* cos(wd.*t);

vud = -gamma.*exp(-gamma.*t).*A.*cos(wd.*t)-exp(-gamma.*t).*A.*sin(wd.*t).*wd;
E = (1/2).*m.*vud.^2 + (1/2).*k.*xud.^2;
dEdt = -c.*vud.^2;
plot(t,xud,'b');
plot(t,vud,'r');
plot(t,E./E(1),'g');
plot(t,dEdt,'c');

m=m1*10;k=k1;c=c1;w0=sqrt(k/m);A=1;gamma=c/2/m;
wd=(w0^2-gamma^2)^0.5;
xud = exp(-gamma.*t) .* A .* cos(wd.*t);
%plot(t,xud,'g');
m=m1;k=k1*10;c=c1;w0=sqrt(k/m);A=1;gamma=c/2/m;
wd=(w0^2-gamma^2)^0.5;
xud = exp(-gamma.*t) .* A .* cos(wd.*t);
%plot(t,xud,'r');
m=m1;k=k1;c=c1*10;w0=sqrt(k/m);A=1;gamma=c/2/m;
wd=(w0^2-gamma^2)^0.5;
xud = exp(-gamma.*t) .* A .* cos(wd.*t);
%plot(t,xud,'c');
end

%figure(3);
%clf;
%hold on;
%plot(t,xod,'m');
%plot(t,xcd,'r');
%plot(t,xud,'g');
%plot(t,exp(-gamma.*t) .* A,'b');
%plot(t,exp(-gamma.*t) .* -A,'b');
%grid on;

