%% Frequency response of a MEMS Duffing oscillator
clear all; clf;
'//////////////////////////////'
% Given values
%M*x'' + D*x' + K*x + kappa*x = F0*cos(w*t)
%x'' + 2*D_M2w0*w0*x' + w0^2*x + kappa_M*x = F0_M*cos(Omega*t)

if 1
    M = 8e-10; %mass 8e-10
    D = 1.55e-7; %damping 1.55e-7
    K = 2; %stiffness 2
    kappa = 4.1e12; %nonlinearity 4.1e10
    N = 70; e0 = 8.854e-12; h = 10e-6; g = 3e-6; V = 0.9;
    F0 = 2*N*1/2*e0*h/g*V^2; %Force amplitude
    M,D,K,kappa,F0
    w0 = sqrt(K/M)
    Q = w0*M/D
    xmax = F0/w0/D
    amplitude = linspace(0.01, 2, 1000)*xmax; %Range of amplitude
else
    M = 1; %mass
    D = 0.05*2; %damping
    K = 1; %stiffness
    kappa = -1; %nonlinearity
    F0 = 0.2; %Force amplitude
    amplitude = linspace(0.01, 3, 1000); %Range of amplitude
end
w0 = sqrt(K/M); %Primary resonance
kappa_M = kappa/M; 
F0_M = F0/M; 
D_M2w0 = D/(M*2*w0); %linear damping
j1 = 0; j2 = 0; j3 = 0;

%% Finding the values of excitation frequency w/w0= F and G. And corresponding eigen values
for j = 1 : 1 : length(amplitude)
     w_w0_23(j)=(1+3*kappa_M/(8*w0^2)*amplitude(j)^2+sqrt((F0_M/(2*w0^2*amplitude(j)))^2-D_M2w0^2));
     lambda_w_w0_23_1 = sqrt(-(w_w0_23(j)-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_23(j)-w0-9*kappa_M/(8*w0)*amplitude(j)^2))-(D_M2w0*w0);
     %lambda_w_w0_23_2 = -sqrt(-(w_w0_23(j)-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_23(j)-w0-9*kappa_M/(8*w0)*amplitude(j)^2))-(D_M2w0*w0);
     w_w0_1(j)=(1+3*kappa_M/(8*w0^2)*amplitude(j)^2-sqrt((F0_M/(2*w0^2*amplitude(j)))^2-D_M2w0^2));
     lambda_w_w0_1_1 = sqrt(-(w_w0_1(j)-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_1(j)-w0- 9*kappa_M/(8*w0)*amplitude(j)^2))-(D_M2w0*w0);
     %lambda_w_w0_1_2 = -sqrt(-(w_w0_1(j)-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_1(j)-w0- 9*kappa_M/(8*w0)*amplitude(j)^2))-(D_M2w0*w0);
	if lambda_w_w0_23_1==conj(lambda_w_w0_1_1) %Terminate the loop
        break
    else
        if kappa_M > 0 %Spring hardening 
            j1 = j1 + 1;
            X1(j1) = w_w0_1(j);
            Y1(j1) = amplitude(j);
            stability = (w_w0_23(j)*w0-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_23(j)*w0-w0- 9*kappa_M/(8*w0)*amplitude(j)^2)+(D_M2w0*w0)^2;
            if stability < 0
                j2 = j2 + 1;
                X2(j2) = w_w0_23(j);
                Y2(j2) = amplitude(j);
            else
                j3 = j3 + 1;
                X3(j3) = w_w0_23(j);
                Y3(j3) = amplitude(j);
            end
        else %Spring softening
            j1 = j1 + 1;
            X1(j1) = w_w0_23(j);
            Y1(j1) = amplitude(j);
            stability = (w_w0_1(j)*w0-w0-3*kappa_M/(8*w0)*amplitude(j)^2)*(w_w0_1(j)*w0-w0- 9*kappa_M/(8*w0)*amplitude(j)^2)+(D_M2w0*w0)^2;
            if stability < 0
                j2 = j2 + 1;
                X2(j2) = w_w0_1(j);
                Y2(j2) = amplitude(j);
            else
                j3 = j3 + 1;
                X3(j3) = w_w0_1(j);
                Y3(j3) = amplitude(j);
            end
        end
	end
end
%xlim([0,2.5]); % Range of x-axis
%ylim([0,2.5]); %Range of y-axis
figure(1); clf; hold on; grid on;
if exist('X2')
	plot(X2*w0,Y2,'g'); 
end
plot(X1*w0,Y1,'b'); 
plot(X3*w0,Y3,'b');
title(sprintf('M=%0.2g, D=%0.2g, K=%0.2g, kappa=%0.2g, Q=%0.3g',M,D,K,kappa,Q));
ylabel('Displacement amplitude [meters]');
xlabel(strcat('Frequency [rad/sec]', sprintf('      w0=%0.3g, F0=%0.3g',w0,F0)));

 


