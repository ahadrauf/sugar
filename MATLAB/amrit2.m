u = 1e-6; %Micro
lw = 10*u; %linewidth
V = 100; %Applied voltage

H = 10*u; %Thickness
Ndr = 60; %Number of drive fingers
Nsn = 100; %Number of sense fingers 
Nfb = 100; %Number of feedback fingers
Lcomb = 100*u; %Comb finger length

%10um feature size
if lw == 10*u
    L = 1000*u; %Flexures length
    w = lw; g = lw; %Width and gap
    Lm = 11010*u; Wm = 400*u; %Mass geometry
elseif lw == 20*u
    L = 2000*u; %Flexures length
    w = lw; g = lw; %Width and gap
    Lm = 12100*u; Wm = 600*u; %Mass geometry
elseif lw == 40*u
    L = 2000*u; %Flexures length
    w = lw; g = lw; %Width and gap
    Lm = 12390*u; Wm = 600*u; %Mass geometry
end

%% Mass
rho = 2350; %density of Silicon
e0 = 8.85e-12;% Permittivity of Air
E = 169e9; %Elastic modulus
Fdr = Ndr * e0 * H * V^2 / g; %Combdrive force computation
Fsn = Nsn * e0 * H * V^2 / g; %Combdrive force computation
Ffb = Nfb * e0 * H * V^2 / g; %Combdrive force computation
main_mass = rho * H * Lm * Wm;
N = Ndr + Nfb + Nsn;
finger_mass = rho * H * 2*N * w * Lcomb; 
M = main_mass + finger_mass; %Computing the Bulk Mass

%% Stiffness
K = 2 * 192 * E * H * w^3 / (12 * (2*L)^3);  %Four beams combined into two fixed-fixed beams 
omega = sqrt(K/M); %Resonace [rad/sec]

%% Damping
overlap = 40*u; %overlap of fingers
fluid_layer = 10*2*u; %thickness of air under proof mass
mu = 1.75e-5; %viscosity of air [sPa]
Df = 2*N * mu * overlap * H / g; %Damping contribution by the fingers
DM = 2 * mu * Lm * Wm / fluid_layer; %Damping contribution by the proof mass
D = DM + Df; %Total damping
    
%% Amplitude at resonance
'///////////////////////////'
A_dr = Fdr / (D * omega) %Amplitude of Resonance
x_dr = Fdr / K %Static deflection
freq_Hz = omega/2/pi %resonance frequency
Q = omega * M / D %Quality factor computation

%% gravity
Kg = 2 * 192 * E * w * H^3 / (12 * (2*L)^3);  %Four beams combined into two fixed-fixed beams 
x_g = M * 9.8 / Kg
z_freq = sqrt(Kg/M) / 2 / pi

