function de 
'-------------' 
clear all
%DE: Mx'' + Dx' + Kx = Fsin(10 t) 
%IC: x(0)=0, x'(0)=0 
tic
T2 = 0.09e-3;
tspan = 0:T2/10000:T2; %time scale  
initial_x    = 0; 
initial_dxdt = 0; 
 
u = 1e-6; %micron
density = 2300; %density of si [kg/m^3]
length = 100*u; %length of flexure [m]
width = 2*u; %width of flexure [m]
thickness = 2*u; %thickness of flexure [m]
I = thickness * width^3 / 12; %second moment of area
E = 160e9; %Modulus [Pa]
air_damping = 1.5e-7;
K1 = 3*E*I/(length/2)^3; %half flexure
K2 = K1/2; %full flexure
K = K2*2; %2 flexures
M = density * length * width * thickness * 27;
D = air_damping / 1000;
e0 = 8.854e-12; %permittivity of the medium
M,D,K

V = 0.02;
%FF = 1/2 * e0 * length * thickness; %force coefficient, V^2 / gap^2;
FF = 2*25*1/2 * e0 * thickness; %force coefficient, V^2 / gap;
%xmax = (FF * 50^2 / (2*u)^2) / D / sqrt(K/M)
xmax = (FF * V^2 / (2*u)) / D / sqrt(K/M) / 2
tau = 1/(D / 2 / M)
%initial_x    = 3e-6; 

MDKF = [M,D,K,FF,V];
tol = 1e-6;
options = odeset('RelTol',tol,'AbsTol',tol);
[t,x] = ode45( @eqn, tspan, [initial_x initial_dxdt], options, MDKF ); 

figure(2);clf; grid on; hold on;
%plot(t/1e-3,x(:,2)/1e-6); 
plot(t/1e-3,x(:,2)); 
%xlabel('Time [ms]'); ylabel('Displacement [um]'); 
xlabel('Time [ms]'); ylabel('Velocity [m/s]'); 
toc

VDC = 0.09;
VAC = V;
gap = 2*u; 
A = 0;
n = 0;
for j = 2 : size(t,1)-1
    dt = t(j) - t(j-1);
    n = n + 1;
    v(n) = VDC + VAC * sin( 0.5 * omega * t(j) );
    tt(n) = t(j);
    i(n) = v(n) * e0 * length * thickness * x(j,2) / (gap - x(j,1))^3;
    A = A + i(n)^2 * dt;
end
i_rms = sqrt( A / 0.091e-3 )
figure(2);clf; grid on; hold on;
plot(tt/1e-3,i/1e-6); 
xlabel('Time [ms]'); ylabel('Current [uA]'); 
figure(3);clf; grid on; hold on;
plot(tt/1e-3,v); 
xlabel('Time [ms]'); ylabel('Voltage [V]'); 
plot(tt/1e-3,v.^2,'r'); 

K=5.12;gap0=2e-6; h_ox=1e-9; L=100e-6; h=2e-6; e0=8.854e-12; V = sqrt(  K * ( gap0 - h_ox ) / ( 1/2 * e0 * L * h / h_ox^2 ) )


function dxdt = eqn(t,x,MDKF)
    V = MDKF(5);
    K = MDKF(3);
    D = MDKF(2);
    M = MDKF(1);
    FF = MDKF(4); %4.427e-016
    G = 2e-6;
    gap = ( G - x(1) );
    VDC = 0.09;
    %VDC = 0.0895;
    omega = sqrt(K/M);
    
    F = FF * (VDC+V*sin(0.5*omega*t))^2 / gap^2;        
    %F = FF * (V*sin(0.5*omega*t))^2 / (2e-6)^2;
    %F = FF * V^2 / 2e-6 * (sin(0.5*omega*t))^2;
    
    if gap <= 1e-9
        x(1) = 2e-6;
        x(2) = 0;
    end
        
    dxdt_1 = x(2); 
    dxdt_2 = M \ (-D*x(2) - K*x(1) + F); 
    dxdt = [dxdt_1; dxdt_2]; 
end 

end


