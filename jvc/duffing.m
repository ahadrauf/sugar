%figure(1);duffing(.5,.6,1,2) %w0,u,k,aplha
function duffing(w0,u,k,alpha) %w0,u,k,aplha
amax=k/(2*w0*u)
figure(1);
clf;
hold on;
i=0;
for a = 0.1 : amax/500 : amax
   i = i + 1;
   omega(i) = a;
   [sig1(i),sig2(i)] = sigma(alpha,w0,u,k,a);
end

plot(real(sig2),omega,'b');
plot(real(sig1),omega,'b');
grid on;
hold off;

figure(2);
clf;
hold on;
i = 0;
for a = 0.1 : amax/500 : amax
   i = i + 1;
   [sig1(i),sig2(i)] = sigma(alpha,w0,u,k,a);
   gamma1(i) = acos(1/4*a*(3*alpha*a^2-8*sig1(i)*w0)/k);
   gamma2(i) = acos(1/4*a*(3*alpha*a^2-8*sig2(i)*w0)/k);
end
plot(real(sig1),real(gamma1),'b');
plot(real(sig2),real(gamma2),'b');
grid on;
hold off;
