format compact
clear all
clc

%% Constant basic value
u = 1e-6; %micro scale
N_total = 500; %total finger numbers
N_sense = 200; %sense finger numbers


%% Design parameters
V = 60;% voltage can be modified to 10, 20, 30 or so.....
linewidth = 20*u; %Minimum linewidth
L = 2000*u; %length of flexure
%N = 70; %Total number of fingers
w = linewidth; %width of flexure
gf = linewidth; %gap between fingers
wf = linewidth; %width of fingers
Lf = 80*u; %Finger length
h = 10*u; %layer thickness

%% Material properties
rho = 2350; %density of Silicon
e0 = 8.854e-12;% Permittivity of Air
E = 169e9; %Elastic modulus of Si

%% Mass
LM1 = 1e-4; %Length of main mass
LM2 = 2e-2; %Width of main mass
main_mass = rho * h * LM1*LM2; %Computing the proof mass
finger_mass = rho * h * N_total * wf * Lf; %Computing the Finger Mass
M = main_mass + finger_mass; %Computing the Bulk Mass

%% Spring & angular frequency
K = 2*3 * E * h*w^3 / (12 * (L/2)^3);  %Stiffness computation
omega = sqrt(K/M); %Resonace [rad/sec]

%% Damping
overlap = 40*u; %overlap of fingers
h_SiO2 = 10*2*u; %thickness of air under proof mass
mu = 1.75e-5; %viscosity of air [sPa]
Df = (2 *  N_total * mu * overlap * h / gf) + ( N_total * mu * Lf * wf / h_SiO2) ; %Damping contribution by the fingers
DM = mu * LM1*LM2 / h_SiO2; %Damping contribution by the mass
D = DM + Df %Total damping
    
%% Feedback Circuit parameter
gain = 3000; %current amplifier gain
R = 60000; %Resistor to convert current to voltage/ohm.
V_bias = 20; %bias voltage applied on feedback side.



%% Drive force computation:
i = 0;
for N_drive =  70: 1 : 150  %Varying of drive finger numbers
  i=i+1;
%% Drive Force
F_drive =  N_drive * 1/2 * e0 * h * V^2 / gf; %Combdrive force computation

%% Feedback Force
N_feedback = N_total - N_drive - N_sense; %feedback finger numbers
I_sensed = pod_dde_2b_2(N_drive,M,D,h,gf,V,K,N_sense,N_total); %Sensed current.
V_feedback = sqrt(I_sensed * gain)*R+V_bias; %feedback voltage
F_feedback = (N_feedback * 1/2 * e0 * h / gf) * (V_feedback)^2; %feedback force

%% Amplitude at resonance
A = (F_drive-F_feedback) / (D * omega); %Amplitude of Resonance after substracting feedback force. 

%% Damping force calculation
F_damping = D * omega * A; %Damping force

%% Force and finger numbers
x1(i) = N_drive; %Drive finger numbers
x2(i) = N_feedback; %Feedback finger numbers
y1(i) = F_damping; %Damping force
y2(i) = F_feedback; %Feedback force
y3(i) = A; %Amplitude
end



%% Display sensed current, feedback voltage, Amplitude, resonate frequency.
max_sensed_current = max(max(I_sensed)) %Maximum sensed current.
max_feedback_voltage = max(max(V_feedback))-V_bias %Maximum feedback voltage.
max_Amplitude = max(max(A)) %Maximum amplitude
reasonate_frequency = omega %Resonate frequency.




%% Plot damping force & feedback force combined in one window.
figure(1);clf; grid on; plot(x1,y1); %Plot Damping force vs Number of drive Fingers
hold on; plot(x1,y2,'g'); %Plot feedback force vs number of drive fingers 
ylabel('Force[N]'); xlabel('Number of drive Fingers'); title('Damping force & Feedback force'); %title
legend('DampingForce', 'FeedbackForce'); %label

%% Plot damping force vs drive finger numbers individually
%figure(1);clf;hold on; grid on; plot(x1,y1); %Plot Damping force vs Number of drive Fingers
%ylabel('Damping Force [N]'); xlabel('Number of drive Fingers'); title('Damping Force vs Number of Drive Fingers');

%% Plot feedback force vs feedback finger numbers individually
% figure(2);clf; hold on; grid on; plot(x1,y2); %%Feedback Force vs Number of Fingers
% ylabel('Feedback Force [N]'); xlabel('Number of Fingers'); title('Feedback Force vs Number of Feedback Fingers');

%% Plot Amplitude vs drive finger numbers individually
% figure(3);clf; hold on; grid on; plot(x1,y3); %%Amplitude vs Number of Fingers
% ylabel('Amplitude [m]'); xlabel('Number of Fingers'); title('Amplitude vs Number of Drive Fingers');



