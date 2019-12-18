clear all
x = (-1.9: 0.001: -1)*1e-6;
%y = 1./x.^6;
gap = 2e-6;
y = -.1e-6 + 1e-12 * 1e-49 ./abs(gap+x).^7;
figure(3); clf; hold on; grid on;
plot(x,y);
