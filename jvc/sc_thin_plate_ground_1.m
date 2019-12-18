%g=1e-6; L=40e-6; tic;Qthinsc = sc_thin_plate_ground_1( -L + g *i, -0.0001e-6 + g *i);toc
%tic;Qthinsc = sc_thin_plate_ground_1( -10e-6 + 1e-6 *i, -0.001e-6 + 1e-6 *i),toc
%Qe2d=sum(-Q{1}(ceil([min(rhs):max(rhs)+1])))/2,ratio=real(Qthinsc/Qe2d)

function [Q,Cpp,C,Cppcor,Cjvc] = sc_thin_plate_ground_1(za,zb)

%g = 1e-6; %gap between plate and ground
g = abs(imag(zb));
h = 20e-6; %out-of-plane depth
L = abs(abs(real(zb)) - abs(real(za)));
permittivity = 8.854e-12; %permittivity of free space
V = 15/2; %voltage
%za=;
%zb=;

i=0;
for z = [za,zb];
    i=i+1;
%    t(i) = exp(-(g+g*LambertW(exp(-(g-pi*z)/g))-pi*z)/g);
    t(i) = exp(-(1 + LambertW(exp(-(g-pi*z)/g)) - pi*z/g));
end

Q = permittivity * h / pi * V * log(t(2)/t(1));




%D=0e-6;sc_thin_plate_ground_1( -(D+20e-6) + i*2e-6, -D + i*2e-6)

Cpp = permittivity * h * L / (2*g);
C = Q / (2*V);
Cjvc = permittivity * h / pi * log(t(2)/t(1))/2;
Cppcor = permittivity * h * (L+(g/pi)) / (2*g) ;


Cppcor2 = Cpp + permittivity * h / pi / 2;

Qpp = Cpp * (V*2);
Q;
Qppcor = Cppcor * (V*2);

deltaC = C - Cpp;

