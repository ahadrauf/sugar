%corner compliance. x=distance from corner. y=beam deflection w 25uN. beam = 50um2um2um. base = 50um2um100um. 
i=0; i=i+1;
x(i)=50e-6;    y(i)=5.03956e-6; i=i+1;
x(i)=40e-6;    y(i)=5.03958e-6; i=i+1;
x(i)=30e-6;    y(i)=5.03968e-6; i=i+1;
x(i)=20e-6;    y(i)=5.03996e-6; i=i+1;
x(i)=10e-6;    y(i)=5.04103e-6; i=i+1;
x(i)=6e-6;     y(i)=5.04161e-6; i=i+1;
x(i)=4e-6;     y(i)=5.04419e-6; i=i+1;
x(i)=2e-6;     y(i)=5.05423e-6; i=i+1;
x(i)=1e-6;     y(i)=5.07929e-6; i=i+1;
x(i)=0e-6;     y(i)=5.19375e-6; i=i+1;

%base compliance. a=width of base to beam. y=beam deflection w 25uN. middle beam = 50um2um2um. base = 50um2umXum. 
i=0; i=i+1;
a(i)=50e-6;    b(i)=5.03956e-6; i=i+1;
a(i)=10e-6;    b(i)=5.03651e-6; i=i+1;
a(i)=8e-6;     b(i)=5.03557e-6; i=i+1;
a(i)=6e-6;     b(i)=5.03366e-6; i=i+1;
a(i)=4e-6;     b(i)=5.02874e-6; i=i+1;
a(i)=2e-6;     b(i)=5.00713e-6; i=i+1;
a(i)=1e-6;     b(i)=4.98439e-6; i=i+1;
a(i)=0e-6;     b(i)=4.86809e-6; i=i+1;
                      
%plot of corner compliance
figure(1);clf;
plot(x(1:10),y(1:10),'-*');
xlabel('Canilever distance from corner, [m].      Simulated on IntelliSuite. Max automesh=2um');
ylabel('Max deflection with a force of 25uN, [m].  ');
title('Corner Compliance. beam(xyz)=50um2um2um. base(xyz)=50um100um2um. F=25uN. E=160Pa. Nu=0.226');
grid on;

%plot of corner compliance convergence
figure(2);clf;
plot(x(1:9),-y(1:9)+y(2:10),'-*');
xlabel('Cantilever distance from corner, [m].      Simulated on IntelliSuite. Max automesh=2um');
ylabel('y2-y1 deflection with a force of 25uN, [m]');
title('Corner Compliance Convergence. beam(xyz)=50um2um2um. base(xyz)=50um100um2um. F=25uN. E=160Pa. Nu=0.226');
grid on;

%plot of base compliance
figure(3);clf;
plot(a(1:8),b(1:8),'-*');
xlabel('Width of base, [m].      Simulated on IntelliSuite. Max automesh=2um');
ylabel('Max deflection with a force of 25uN, [m]');
title('Base Compliance. beam(xyz)=50um2um2um. base(xyz)=50um100um2um. F=25uN. E=160Pa. Nu=0.226');
grid on;

y1um=4.87864e-6
y2um=4.86809e-6
F=25e-6;
E=160e9;
h=2e-6;w=2e-6;
I=w^3*h/12;
L=50e-6;
Roark=F*L^3/3/E/I %4.8828125e-006
error1um=(Roark-y1um)/Roark
error2um=(Roark-y2um)/Roark

