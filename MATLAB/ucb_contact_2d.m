%Constants
    clear all; format compact; format long g;
    K = 0.0941780821917766; %Stiffness
    gap = 2e-6; %Initial gap between plates
    permittivity = 8.854e-12; %Permittivity of the medium
    length = 100e-6; %Overlapping length of the plates
    thickness = 20e-6; %layer thickness of the plates
    B = permittivity*length*thickness; %Constant
    V = 3.56; %Applied voltage difference between plates
    tic
    
%Contact parameter
    x0 = 1.9e-6;
    %U_spring = (1/2*K*x0.^2);
    %U_capacitor = (-1/2*B*V^2 ./ (gap-x0));
    %C = (U_spring + U_capacitor) / (gap-x0).^-6;
    
%Force
    x = [-0.1:0.0001:1.95]*1e-6; 
    C = (K*x0 - (1/2*B*V^2)/(gap-x0)^2) / (6*(gap-x0)^-7);
    F_spring = -K*x;
    F_capacitor = 1/2*B*V^2./(gap-x).^2;
    F_contact = 6*C*abs(gap-x).^-7;
    figure(2); 
    %clf; 
    grid on; hold on;
    %plot(x, F_contact, 'b');
    %plot(x, F_spring, 'r');
    %plot(x, F_capacitor, 'm');
    plot(x, F_spring +  F_capacitor + F_contact, 'r'); title('F(x)'); xlabel('x'); ylabel('F');
    
%Potential
    U_spring = 1/2*K*x.^2;
    U_capacitor = -1/2*B*V^2./(gap-x);
    U_contact = -C*(gap-x).^-6;
    figure(10); 
    %clf; 
    grid on; hold on;
    %plot(x, U_contact, 'b'); plot(x, U_spring, 'r'); plot(x, U_capacitor, 'm'); 
    plot(x, U_spring +  U_capacitor + U_contact, 'r'); title('U(x)'); xlabel('x'); ylabel('U');

j = 0;
for V = 0 : .001: 5    
    j = j + 1; %plot data index
    

%Newton solver
    C = (-K*x0 - (1/2*B*V^2)/(gap-x0)^2) / -(6*(gap-x0)^-7);
    x = gap * 0; %Initial guess
    delta_x = 0.00001e-6;
    iter = 0;
    while (abs(delta_x)) >= 1e-12 && iter < 200
      x = x - delta_x;
      F =  -K*x -    1/2*B*V^2*(gap+x)^-2 -     (6*-C*(gap+x).^-7 ); %sum of forces to zero out to find x
      dF = -K   - 2*1/2*B*V^2*(gap+x)^-3 -  -7*(6*-C*(gap+x).^-8 ); %dF = slope = delta_F / delta_x
      %F = -K*x - 1e-9;  dF = -K;
      delta_x = dF\F;
      iter = iter + 1;
    end
    FF(j) = K*x   +   1*-1/2*B*V^2./(gap-x).^2  +   (-6*-C*(gap-x).^-7);
%    x, V
    VV(j) = V;
    xx(j) = x;
end
figure(4); clf; grid on; hold on; plot(VV,xx); 
title('Displacement vs voltage'); xlabel('Voltage [V]'); ylabel('Displacement [m]');
%figure(5); clf; grid on; hold on; plot(xx,FF); 
%title('x vs F'); xlabel('x'); ylabel('F');
%iter
x
F
toc

