format compact
clear all

i = 0;
%for N = 50 : 1 : 500 
for linewidth = [2 : 0.1 : 10]*1e-6 
i=i+1;

%% Design parameters
V = 100;% voltage
u = 1e-6;
%linewidth = 5*u; %Minimum linewidth
w = linewidth; %width of flexure
gf = linewidth; %gap between fingers
wf = linewidth; %width of fingers
Lf = 60*u; %Finger length
h = 20*u; %layer thickness
N = 70; %Total number of fingers
L = 500*u; %length of flexure

%% Material properties
rho = 2350; %density of Silicon
e0 = 8.85e-12;% Permittivity of Air
E = 169e9; %Elastic modulus

%% Force
F = N * e0 * h * V^2 / gf; %Combdrive force computation
 
%% Mass
LM = 1e-3; %Length of main mass
main_mass = rho * h * LM^2; 
finger_mass = rho * h * 2*N * wf * Lf; 
M = main_mass + finger_mass; %Computing the Bulk Mass

%% Spring & angular frequency
K = 3 * E * h*w^3 / (12 * (L/2)^3);  %Stiffness computation
omega = sqrt(K/M); %Resonace [rad/sec]

%% Damping
overlap = 40*u; %overlap of fingers
h_SiO2 = 2*u; %thickness of air under proof mass
mu = 1.75e-5; %viscosity of air [sPa]
Df = mu * overlap * h / gf; %Damping contribution by the fingers
DM = mu * LM^2 / h_SiO2; %Damping contribution by the mass
D = DM + Df; %Total damping
    
%% Amplitude at resonance
A = F / (D * omega); %Amplitude of Resonance
Q = omega * M / D; %Quality factor computation

y(i) = A;
%x(i) = N;
x(i) = linewidth;

end

figure(1); clf; hold on; grid on; plot(x,y,'-');
ylabel('A');
%xlabel('N');
xlabel('linewidth');

%% Results
'------ Results -------'
F %comb drive force
K %flexure stiffness
M %proof mass
omega %angular resonance frequency
A %amplitude of motion

