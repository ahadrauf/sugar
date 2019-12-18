function integrand=sc9(t,g)
integrand = 1./t.*pi./(g.*(1.+1./t));

H = 10e-6;
L = 20e-6;
V = 50;
A = H * L;
E = 8.854e-12;

integrand = V * H * E / 2/pi * integrand;
Fpp = 0.5 * E * A * V^2 / g^2;

