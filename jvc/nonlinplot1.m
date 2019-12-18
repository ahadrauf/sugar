%nonlinplot1

x=0:0.5/1000:0.4;

ylin=3*x;
ynonlin=3*x+3.26*x.^3;
figure(1);
clf;
plot(ylin,x,'r',ynonlin,x,'b');
grid on;
