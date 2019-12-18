
%Constants
    clear all; format compact; format long g;
    K = 0.1; %Stiffness
    gap = 2e-6; %Initial gap between plates
    permittivity = 8.854e-12; %Permittivity of the medium
    length = 2e-6; %Overlapping length of the plates
    thickness = 20e-6; %layer thickness of the plates
    B = permittivity*length*thickness; %Constant
    V = 30; %Applied voltage difference betwee plates
        
%Contact parameter
    x0 = 1.96e-6;
    U_spring = (1/2*K*x0.^2);
    U_capacitor = (-1/2*B*V^2 ./ (gap-x0));
    C = (U_spring + U_capacitor) / (gap-x0).^-6;

%Potential
    x = [-1:0.0001:1.955]*1e-6; 
    U_spring = 1/2*K*x.^2;
    U_capacitor = -1/2*B*V^2./(gap-x);
    U_contact = -C*(gap-x).^-6;
    figure(1); clf; grid on; hold on;
    plot(x, U_contact, 'b'); plot(x, U_spring, 'r'); plot(x, U_capacitor, 'g'); 
    plot(x, U_spring +  U_capacitor + U_contact, 'k'); xlabel('Displacement x [m]'); ylabel('Potential energy U [J]');
    title('Uspring(red)    Ucontact(blue)   Ucapacitor(green)   Unet(black)'); 
%Force
    x = [-20:0.0001:1.95]*1e-6; 
    F_spring = -K*x;
    F_capacitor = 1/2*B*V^2./(gap-x).^2;
    F_contact = 6*C*abs(x-2e-6).^-7;
    figure(2); clf; grid on; hold on;
    plot(x, F_contact, 'b');plot(x, F_spring, 'r');plot(x, F_capacitor, 'g');
    plot(x, F_spring +  F_capacitor + F_contact, 'k'); 
    title('Fspring(red)    Fcontact(blue)   Fcapacitor(green)   Fnet(black)'); 
    xlabel('x'); ylabel('F');

%Newton solver
    x = gap / 3; %Initial guess
    delta_x = -0.00001e-6;
    iter = 0;
    while (abs(delta_x)) >= 1e-20 && iter < 1000
      x = x - delta_x;
      F =  -K*x +    -1*-1/2*B*V^2./(gap-x).^2 +    (-6*-C*(gap-x).^-7);
      dF = -K   + -2*-1*-1/2*B*V^2./(gap-x).^3 + -7*(-6*-C*(gap-x).^-8);
      delta_x = dF\F;
      iter = iter + 1;
    end
    
    disp('---------- Results ----------')
    delta_x
    iter
    dF
    F
    x
