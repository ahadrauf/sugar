function xdot=vacdamp1(t,x)

m=1;
c=1;
k=1;
w0=sqrt(k/m);
F0=1;

F=F0*cos(w0*t) + 0* c*x(2);


xdot(1) = x(2);
xdot(2) = 1/m * (F - c*xdot(1) - k*x(1));

xdot = [x(2);xdot(2)];

