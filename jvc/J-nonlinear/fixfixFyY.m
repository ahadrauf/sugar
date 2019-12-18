function [Fy,Y,M]=fixfixFy(Fx,k)
h=6e-6;
w=6e-6;
L=100e-6;
E = 185e9;
I = h * w^3 / 12;
A = h * w;
lam = L / 2 * sqrt( Fx / E / I );

N1 = lam^3 * 8 * sqrt(2);
D = sqrt(3/2 - 1/2 * (tanh(lam))^2 - 3/2 * tanh(lam) / lam);

Fy = sqrt(I/A) * E * I / L^3  * N1 / D;

N2 = sqrt(I/A) * 2 * sqrt(2) * (lam - tanh(lam));
Y = N2/D;

N3 = sqrt(I/A) * E*I / L^2 * lam^3 * 4 * sqrt(2) * tanh(lam);
M = N3/D;

%[f,y]=fixfixFyY(300e-4)
%f = 5.223913491604605e-003
%y = 8.790989065775853e-006

%[f,y]=fixfixFyY(100e-4)
%f = 1.806269814656536e-003
%y = 5.030236143960421e-006

%[f,y]=fixfixFyY(50e-4)
%f = 1.062770294762911e-003
%y = 3.547462259020234e-006
%m = 3.501590893889557e-008

%k = 0.69;
b = 12 * E * I / L^3 / k;
c = -Fy / k;
x = (-c/2 + sqrt(c^2/4+b^3/27))^(1/3) - (c/2 + sqrt(c^2/4+b^3/27))^(1/3);
x

