function harm_phase

w0=1;
m=1;
c=(1/2)*m*w0;
gamma=c/2/m;
figure(1);%clf;
i=0;
for w=0:0.01:3
i=i+1;
phase(i)=atan(2.*gamma.*w./(w0.^2-w.^2));
if phase(i)<0
    phase(i)=phase(i)+pi;
end
W(i)=w;
end
plot(W,phase);
grid on;
hold on;