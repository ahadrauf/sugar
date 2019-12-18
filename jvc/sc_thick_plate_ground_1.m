%D=0.001e-6;DD=1e-6;Qthin=sc_thin_plate_ground_1( -(D+DD*1e-6) + i*2e-6, -D + i*2e-6),Qthick=sc_thick_plate_ground_1( -(D+DD*1e-6) + i*2e-6, -D + i*2e-6),percent=real(Qthick/Qthin*100)

% tic;Qthicksc = sc_thick_plate_ground_1( -10e-6 + 1e-6 *i, -0.001e-6 + 1e-6 *i),toc,Qe2d=sum(-Q{1}(ceil([min(rhs):max(rhs)+1])))/2,ratio=real(Qthicksc/Qe2d)

function Q = sc_thick_plate_ground_1(za,zb)

G = 1e-6; %gap between plate and ground
P = 1e-6; %plate thickness
h = 20e-6; %out-of-plane depth
permittivity = 8.854e-12; %permittivity of free space
V = 15/2; %voltage
%za=;
%zb=;

i=0;
g = strcat('g = ',num2str(G));
p = strcat('p = ',num2str(P));
for z = [za,zb];
    Zr = real(z);
    Zi = imag(z);    
    Z = strcat('Z = ',num2str(Zr),' + ',num2str(Zi), ' * sqrt(-1) ' );
    
    i=i+1;
    [t]=solve('Z =g*((2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1))*arctanh(sqrt((t+1)/(t-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1))))/(sqrt(-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1)))+(-2+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1))*sqrt((t+1)/(t-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1)))/(sqrt(-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1))*(1-(t+1)/(t-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1))))+ln((sqrt((t+1)/(t-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1)))*sqrt(-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1))-1)/(sqrt((t+1)/(t-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1)))*sqrt(-1+2*(1+p/g)^2+2*(1+p/g)*sqrt((1+p/g)^2-1))+1)))/3.1415926535897',p,g,Z);
    T(i)=t.t;
end

Q = permittivity * h / pi * V * log(T(2)/T(1));

%D=0e-6;sc_thick_plate_ground_1( -(D+20e-6) + i*2e-6, -D + i*2e-6)

Cpp = permittivity * h^2 / G;
Qpp = Cpp * V;

