%% Frequency response of a MEMS Duffing oscillator
clear all; clf; 

%Equations:
%M*x'' + D*x' + K*x + kappa*x = F0*cos(w*t)
%x'' + 2*D_M2w0*w0*x' + w0^2*x + kappa_M*x = F0_M*cos(Omega*t)
'/////////////////////////////'
E=160e9; 
L=294.7e-6; H=20e-6; W=2e-6;
I=W^3*H/12;
alpha = 2.5e-6;
DT = 100;
K1 = 3*E*I/(L/2)^3
L=L*(1+alpha*DT);
H=H*(1+alpha*DT);
W=W*(1+alpha*DT);
I=W^3*H/12;
K2 = 3*E*I/(L/2)^3
M = 8e-10; 
w1=sqrt(K1/M)
w2=sqrt(K2/M)
K=K2


%Parameters
    M = 8e-10; %mass 8e-10
    D = 1.55e-7; %damping 1.55e-7
    K = 2; %stiffness 2
    kappa = -1000 * 4.1e12; %nonlinearity 4.1e10
    N = 70; e0 = 8.854e-12; h = 10e-6; g = 3e-6; V = 0.9;
    F0 = 2*N*1/2*e0*h/g*V^2; %Force amplitude
    w0 = sqrt(K/M);
    Q = w0*M/D;
    xmax = F0/w0/D;
    amplitude = linspace(0.01, 2, 1000)*xmax; %Range for displacement amplitude
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
    figure(1); hold on; grid on;
    if exist('X2') %If there is an unstable curve
        plot(X2*w0,Y2,'g'); 
    end
    plot(X1*w0,Y1,'b'); 
    plot(X3*w0,Y3,'b');
    title(sprintf('M=%0.2g, D=%0.2g, K=%0.2g, kappa=%0.2g, Q=%0.3g',M,D,K,kappa,Q));
    ylabel('Displacement amplitude [meters]');
    xlabel(strcat('Frequency [rad/sec]', sprintf('      w0=%0.3g, F0=%0.3g',w0,F0)));
 
