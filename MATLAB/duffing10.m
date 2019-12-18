%% Frequency response of a MEMS Duffing oscillator
clear all; 
%clf; 
tic
'////////////////////////'
%Equations:
%M*x'' + D*x' + K*x + kappa*x = F0*cos(w*t)
%x'' + 2*D_M2w0*w0*x' + w0^2*x + kappa_M*x = F0_M*cos(Omega*t)
u = 1e-6;
E=160e9; 
L=294.7e-6; H=20e-6; W=2e-6;
I=W^3*H/12;
alpha = 2.6e-6;
DT = 0;
K1 = 3*E*I/(L/2)^3;
L=L*(1+alpha*DT);
H=H*(1+alpha*DT);
W=W*(1+alpha*DT);
I=W^3*H/12;
K2 = 3*E*I/(L/2)^3;
M = 8e-10; 
w1=sqrt(K1/M);
w2=sqrt(K2/M);


%Parameters
    M = 8e-10; %mass 8e-10
    %D = 1.55e-7 %Air
%kappa = 4.2^2 * 37 * 8.854e-12 * 20e-6 / (2e-6 * (5.28e-6)^3)
kappa = 20^2 * 37 * 8.854e-12 * 20e-6 / (2e-6 * (5.28e-6)^3)
    N = 37; e0 = 8.854e-12; 
    V = 1
    g = 2e-6; 
%Temperature    
    DT = 0;
    T = 300 + DT;
    alpha = 2.6e-6;
    E = exp( 2.61 * 1e-3 * 1.6e-19 / T / 1.38e-23) * 145e9
    L=L*(1+alpha*DT);
    H=H*(1+alpha*DT);
    W=W*(1+alpha*DT);
    I=W^3*H/12;
%K = 3*E*I/(L/2)^3 * 0.9
K = 3*E*I/(L/2)^3
    
    w0 = sqrt(K/M)
    %Q = w0*M/D
    Q = 5e4;
    D = w0*M/Q
    %xmax = F0/w0/D
    F0 = 0.5* 2*N*1/2*e0*H/g*V^2; %Force amplitude
    xmax = F0/w0/D
    amplitude = linspace(0.0001, 2, 10000)*xmax; %Range for displacement amplitude
    
if 1
        xmax = 45e-6;
        xdot2 = xmax * w0; 
        w_hour = 2*pi / (60) /2; %Radians per sec
        F0 = 2 * M * w_hour * xdot2; %Coriolis amplitude    
        amplitude = linspace(0.001, 1, 10000)*8e-6; %Range for displacement amplitude
    end
    
    %amplitude = linspace(0.00001, 0.0005, 1000)*xmax; %Range for displacement amplitude
%Eqn conversion 
    w0 = sqrt(K/M); %Primary resonance
    kappa_M = kappa/M;
    F0_M = F0/M; 
    D_M2w0 = D/(M*2*w0); %linear damping
j1 = 0; j2 = 0; j3 = 0;

%% Plots 
for j = 1 : 1 : length(amplitude)
    %Eigen values and w/w0
        w_w0_23(j)=(1+3*kappa_M/(8*w0^2)*amplitude(j)^2+sqrt((F0_M/(2*w0^2*amplitude(j)))^2-D_M2w0^2));
        lambda_w_w0_23_1 = sqrt(-(w_w0_23(j)-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_23(j)-w0-9*kappa_M/(8*w0)*amplitude(j)^2))-(D_M2w0*w0);
        w_w0_1(j)=(1+3*kappa_M/(8*w0^2)*amplitude(j)^2-sqrt((F0_M/(2*w0^2*amplitude(j)))^2-D_M2w0^2));
        lambda_w_w0_1_1 = sqrt(-(w_w0_1(j)-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_1(j)-w0- 9*kappa_M/(8*w0)*amplitude(j)^2))-(D_M2w0*w0);
	if lambda_w_w0_23_1==conj(lambda_w_w0_1_1) 
        break %Terminate 
    else
        if kappa_M > 0 %Spring hardening case
            j1 = j1 + 1;
            X1(j1) = w_w0_1(j);
            Y1(j1) = amplitude(j);
            stability = (w_w0_23(j)*w0-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_23(j)*w0-w0- 9*kappa_M/(8*w0)*amplitude(j)^2)+(D_M2w0*w0)^2;
            if stability < 0 %if < 0, then it is unstable
                j2 = j2 + 1;
                X2(j2) = w_w0_23(j);
                Y2(j2) = amplitude(j);
            else %if stable
                j3 = j3 + 1;
                X3(j3) = w_w0_23(j);
                Y3(j3) = amplitude(j);
            end
        else %Spring softening case
            j1 = j1 + 1;
            X1(j1) = w_w0_23(j);
            Y1(j1) = amplitude(j);
            stability = (w_w0_1(j)*w0-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_1(j)*w0-w0- 9*kappa_M/(8*w0)*amplitude(j)^2)+(D_M2w0*w0)^2;
            if stability < 0 %if < 0, then it is unstable
                j2 = j2 + 1;
                X2(j2) = w_w0_1(j);
                Y2(j2) = amplitude(j);
            else %if stable
                j3 = j3 + 1;
                X3(j3) = w_w0_1(j);
                Y3(j3) = amplitude(j);
            end
        end
	end
end
%Figure displacy
    k = 1000;   
    figure(1); 
%clf;
    hold on; grid on;
    if exist('X2') %If there is an unstable curve
        plot(X2*w0/k/2/pi,Y2,'g'); 
    end
%color = 'b';    
color = 'r';    
    plot(X1*w0/k/2/pi,Y1,color); 
    plot(X3(1:length(X3)-1)*w0/k/2/pi,Y3(1:length(X3)-1),color);
    title(sprintf('PODMEMS control with increased nonlinearity',M,D,K,kappa,Q));
    ylabel('Displacement amplitude [m]');
    %xlabel(strcat('Frequency [kHz]', sprintf('      w0=%0.3g, F0=%0.3g',w0,F0)));
    xlabel('Frequency [kHz]');
toc    
