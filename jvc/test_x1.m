function x1
p=2;
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

f(i)=12*q(i);

end

%plot(x,Fy,'r',x,q,'g',x,M,'b');
plot(Fy,q,f,q,M,q);
