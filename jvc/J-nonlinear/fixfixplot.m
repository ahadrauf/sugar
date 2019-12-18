%fixfixplot(2)
function fixfixplot(p,n)
%p=2;
i = 0;
for lam = 0.001 : p/50 : p
   
i = i + 1;
x(i) = lam;
FyN = lam^3 * 8 * sqrt(2);
FyD = sqrt( 3/2 - (1/2)*(tanh(lam))^2 - 3/2 * (tanh(lam))/lam);
Fy(i) = FyN / FyD;

qN = 2 * sqrt(2) * (lam - tanh(lam));
qD = sqrt( 3/2 - (1/2)*(tanh(lam))^2 - 3/2 * (tanh(lam))/lam);

q(i) = qN/qD;

MN = lam^3 * 4 * sqrt(2) * tanh(lam);
MD = sqrt( 3/2 - (1/2)*(tanh(lam))^2 - 3/2 * (tanh(lam))/lam);
M(i) = MN/MD;

%stiffness for force
kfN = 4 * lam^3;
kfD = lam - tanh(lam);
kf(i) = kfN/kfD;

%stiffness for moment
kmN = 2 * lam^3 * tanh(lam);
kmD = lam - tanh(lam);
km(i) = kmN/kmD;

fl(i)=12*q(i);
ml(i)=6*q(i);

Fx(i) = 4 * lam^2;
Fxl(i) = 0;

end

clf;
hold on;

plot(q,Fy,'r',q,fl,'m',q,M,'b',q,kf,'k',q,km,'g',q,ml,'c', q, Fx, 'b-x',q, Fxl, 'b-x');

%thick plot
%for s=-0.5 : 0.25 : 0.5
%   plot(q,Fy+s,'r',q,fl+s,'m',q,M+s,'b',q,kf+s,'k',q,km+s,'g',q,ml+s,'c');
%end

%TEXT
grid on;
xlabel('q  transverse displacement   [nondimensional]');ylabel('F  Force,  M  Moment,  K  Stiffness   [nondimensional]');
a=0;b=-Fy(i)/20;c=16;
a=a+b;
text(q(i)/c,Fy(i)+a,'F-nonlinear = RED');a=a+b;
text(q(i)/c,Fy(i)+a,'F-linear = MAGENTA');a=a+b;
text(q(i)/c,Fy(i)+a,'M-nonlinear = BLUE');a=a+b;
text(q(i)/c,Fy(i)+a,'M-linear = CYAN');a=a+b;
text(q(i)/c,Fy(i)+a,'K-nonlinear wrt F = BLACK');a=a+b;
text(q(i)/c,Fy(i)+a,'K-nonlinear wrt M = GREEN');
title('Modeling of the transverse deflection of a fixed-fixed beam. Theory');

F=Fy(i)
Q=q(i)
k=(F-12*Q)/Q^3
F=12*Q+k*Q^3
k=0.69
plot(q(1:i),12*q(1:i)+k*q(1:i).^3,'ko');

%n=0.3;
F=M(i)
Q=q(i)
k=(F-6*Q-n*Q^.5)/Q^3
F=6*Q+k*Q^3
k=0.37
plot(q(1:i),6*q(1:i)+k*q(1:i).^3+n*q(1:i).^.5,'ko');

F=Fx(i)
Q=q(i)
k=(F)/Q^2
F=Q^3
k=0.6
plot(q(1:i),k*q(1:i).^2,'ko');

hold off;
