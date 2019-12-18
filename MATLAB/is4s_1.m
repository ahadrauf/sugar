'///////////////// IS4S IMU /////////////////'
clear all

%% Proof mass 
    density = 2300; %kg/m^3, density of silicon
    u = 1e-6; %micro 
    length = 500*u; %m, length of proof mass 
    thickness = 8*u; %m, thickness of proof mass 
    area = length * length; %m^2, area of proof mass
    volume = area * thickness; %m^3, volume of proof mass
    M_proofmass = density * volume %kg, mass of the IMU proof mass
    air_viscosity = 1.78e-5; %air viscosity
    fluid_gap = 2*u; %m, viscous fluid of air under proof mass
    D_proofmass = air_viscosity * area / fluid_gap %N*s^2/m
 
%% Piezo launch stiffness 
    G = 9.8; %m/s^2, acceleration due to gravity
    accel_launch = 50e3*G; %m/s^2, acceleration for launch 50,000*G
    x_piezo_launch = 50*u; %Max amplitude
    K_piezo_launch = 10*M_proofmass * accel_launch / x_piezo_launch %N/m
    
%% Piezo drive force
    F_piezo_drive = K_piezo_launch * x_piezo_launch %N, force required for deflection
    V_piezo_drive = 40; %V, piezo drive voltage
    piezo_V2x = x_piezo_launch / V_piezo_drive; %proportionality btwn x and V

%% Piezo resonance
    w0_piezo = sqrt(K_piezo_launch/M_proofmass) %rad/s, resoance of piezo drive mode
    f0_piezo_drive = w0_piezo / (2*pi); %Hz, resoance of piezo drive mode

%% SOI stiffness
    N = 250; %number of comb fingers 
    eps = 8.854e-12; %permittivity of the medium 
    V = 125; %V, comb drive voltage
    g = 2*u; %m, gap btwn comb fingers
    h = thickness; %m, layer thickness of comb drive
    roll_diameter = 155e-3; %m, diameter of roll 
    r = roll_diameter / 2; %m
    K_soi = 38; %N/m, stiffness of Turner gyro
    w0_soi = sqrt(K_soi / M_proofmass); %rad/s
    F_soi = 2 * 1/2 * N * eps * h * V^2 / g %N, drives primary mode
    x_soi = F_soi / K_soi
    dxdt = w0_soi * x_soi;
    w = 36e3 * pi/180; %rad/s, max roll rate is 36k deg/s
    F_Coriolis = 2 * M_proofmass * w * dxdt %N, Coriolis force, drives secondary mode, cross(w,dxdt)
    F_centrifugal = M_proofmass * w * w * r %N, yaw z-axis (turn left|right), pitch y-axis (slope up|down), roll x-axis (twist), cross(w,cross(w,r))
    ddrdtdt = accel_launch;
    F_transnational = M_proofmass * ddrdtdt %N, launch acceleration
    dwdt = 0; %rad/s^2, changing agular rate
    F_Euler = M_proofmass * dwdt * r %N, for w changing in speed or direction, cross(dwdt,r)
    K_soi = F_centrifugal / x_soi

%% (A) Static displacement vs acceleration launch +/-50k*G  
    a = [-accel_launch : (accel_launch- -accel_launch)/10 : accel_launch];
    x = (M_proofmass * a) / K_piezo_launch;
    figure(1); clf; grid on; hold on;
    plot(a/G/1000,x/u); xlabel('Launch acceleration [k*G*m/s^2]'); ylabel('Deflection of proof mass [um]'); title('Piezo-IMU stiff          Launch accel = +/-50k*G');    
    coefficient_of_friction_of_Si = (0.7-0.2)/2; %
    gap_spacing = 0.1*u;  
    V_hold_down_launch = sqrt( M_proofmass * accel_launch / (eps * area / gap_spacing^2) )
	F_hold_down_launch = eps * area / gap_spacing^2 * V_hold_down_launch^2 / coefficient_of_friction_of_Si
    if (M_proofmass * a) > coefficient_of_friction_of_Si * F_hold_down_launch
        x = (M_proofmass * a) / K_soi;
    else
        x = 0*a; %m, Hold down displacement
    end
    figure(2); clf; grid on; hold on;
    plot(a/G/1000,x/u); xlabel('Launch acceleration [k*G*m/s^2]'); ylabel('Deflection of proof mass [um]'); title('SOI-IMU Lock down          Launch accel = +/-50k*G');    
    
%% (B) Static displacement vs dynamic range accel +/-150*G 
        accel_dynamic = 150 * G; %tranlational acceleration
        a = [-accel_dynamic : (accel_dynamic - -accel_dynamic)/10 : accel_dynamic ];
        x = (M_proofmass * a) / K_piezo_launch;
        figure(3); clf; grid on; hold on;
        plot(a/G,x/u); xlabel('Dynamics acceleration [G*m/s^2]'); ylabel('Deflection of proof mass [um]'); title('Piezo-IMU           Dynamic accel = +/-150*G');    
        a = [-accel_dynamic : (accel_dynamic - -accel_dynamic)/10 : accel_dynamic ];
        x = (M_proofmass * a) / K_soi;
        figure(4); clf; grid on; hold on;
        plot(a/G,x/u); xlabel('Dynamics acceleration [G*m/s^2]'); ylabel('Deflection of proof mass [um]'); title('SIO-IMU           Dynamic accel = +/-150*G');    

%% (E) Temperature -55C to 140C

%% (F) Bandwidth > 70Hz
    w0 = sqrt(K_piezo_launch / M_proofmass)
    M = M_proofmass; 
    D = D_proofmass / 1000;
    i = 0;
	K1 = K_piezo_launch; 
    K2 = 40; 
    w_in = 2*pi/(24*60*60)
    %for w = [9.316 : (-9.316 + 9.332)/10000 : 9.332]*1e4 
    for w = 0 : sqrt(K1/M) * 0.6 / 10000 : sqrt(K1/M) * 0.6;
        i = i + 1;
        F1 = F_piezo_drive;
        x1(i) = F1 / sqrt( (K1 - M*w^2)^2 + (D*w).^2 );
        dxdt = w*x1(i);
        F2 = 2 * M * w_in * dxdt;
        x2(i) = F2 / sqrt( (K2 - M*w^2)^2 + (D*w)^2 );
        W(i) = w;
    end
    figure(7); clf; grid on; hold on;
    k=1e3; hz=1/(2*pi);
    plot(W/k*hz,x1/u,'-b'); xlabel('Drive frequency [kHz]'); ylabel('Amplitude [um]'); title('Sdrive-IMU.   Input disturbance: Angular rate = Earth rotation');    
    plot(W/k*hz,x2/u,'-r');
    w10 = sqrt(K1/M)
    w20 = sqrt(K2/M)
    w = w20;
    x1 = F_piezo_drive / sqrt( (K1 - M*w^2)^2 + (D*w)^2 )
    dxdt1 = w * x1;
    F_Coriolis = 2 * M * w_in * dxdt1
    x2_res = F_Coriolis / (D * w20)
    Q_piezo = w10 * M / D
    BW_piezo_Hz = D / M * 1/(2*pi)

    
%    Q1 = w10 * M / D
%    Q2 = w20 * M / D
%    BW1_piezo_Hz = D / M * 1/(2*pi)
%    BW2_piezo_Hz = D / M * 1/(2*pi)

    
%    w0 = sqrt(K_soi / M_proofmass)
%    K = K_soi; M = M_proofmass; D = D_proofmass;
%    w = [0 : sqrt(K/M)/1000000*2 : sqrt(K/M)*2];
%    V_drive = 15
%    F0 = 2 * 1/2 * N * eps * h / g * V_drive^2
%    x = F0 ./ sqrt( (K - M.*w.^2).^2 + (D.*w).^2 );
%	figure(7); %clf; grid on; hold on;
%    plot(w,x,'-r'); xlabel('angular rate [rad/s]'); ylabel('Amplitude [m]'); title('Piezo/SOI-IMU');    
    
    K = K_soi;
    M = M_proofmass; 
    D = D_proofmass;
    w0 = sqrt(K / M)
    i = 0;
	K1 = K; 
    K2 = K; 
    w_in = 2*pi/(24*60*60)
    %for w = [9.316 : (-9.316 + 9.332)/10000 : 9.332]*1e4 
    for w = 0 : sqrt(K1/M) * 1.5 / 10000 : sqrt(K1/M) * 1.5;
        i = i + 1;
        V_drive = 15;
        F1 = 2 * 1/2 * N * eps * h / g * V_drive^2;
        x1(i) = F1 / sqrt( (K1 - M*w^2)^2 + (D*w).^2 );
        dxdt = w*x1(i);
        F2 = 2 * M * w_in * dxdt;
        x2(i) = F2 / sqrt( (K2 - M*w^2)^2 + (D*w)^2 );
        W(i) = w;
    end
    figure(8); clf; grid on; hold on;
    k=1e3; hz=1/(2*pi);
    plot(W/k*hz,x1/u,'-b'); xlabel('Drive frequency [kHz]'); ylabel('Amplitude [um]'); title('Comb-IMU.   Input disturbance: Angular rate = Earth rotation');    
    plot(W/k*hz,x2/u,'-r');
    w10 = sqrt(K1/M)
    w20 = sqrt(K2/M)
    w = w20;
    x1 = F1 / sqrt( (K1 - M*w^2)^2 + (D*w)^2 )
    dxdt1 = w * x1;
    F_Coriolis = 2 * M * w_in * dxdt1
    x2_res = F_Coriolis / (D * w20)
    x_max = F1 / (D * w10)
    Q_soi = w10 * M / D
    BW_soi_Hz = D / M * 1/(2*pi)
    
%% freq shift
    K = K_soi;
    M = M_proofmass; 
    D = D_proofmass/100;
    w0 = sqrt(K / M)
    i = 0;
	K1 = K; 
    K2 = K; 
    w_in = 2*pi/(24*60*60)
    clear x1 x2 W
    %for w = [9.316 : (-9.316 + 9.332)/10000 : 9.332]*1e4 
    for f = 14585 :0.001: 14595
        %: sqrt(K1/M) * 1.5 / 10000 : sqrt(K1/M) * 1.5;
        w = f * 2 * pi;
        i = i + 1;
        V_drive = 1;
        F1 = 2 * 1/2 * N * eps * h / g * V_drive^2;
        x1(i) = F1 / sqrt( (K1 - M*w^2)^2 + (D*w).^2 );
        dxdt = w*x1(i);
        F2 = 2 * M * w_in * dxdt;
        x2(i) = F2 / sqrt( (K2 - M*w^2)^2 + (D*w)^2 );
        W(i) = w;
    end
    figure(11); clf; grid on; hold on;
    k=1e3; hz=1/(2*pi);
    plot(W/k*hz,x1/u,'-b'); xlabel('Drive frequency [kHz]'); ylabel('Amplitude [um]'); title('Comb-IMU.   Input disturbance: High Q');    
    plot(W/k*hz + 4/k,x1/u,'-r'); 
    plot(W/k*hz,x2/u,'-g');
    w10 = sqrt(K1/M)
    w20 = sqrt(K2/M)
    w = w20;
    x1 = F1 / sqrt( (K1 - M*w^2)^2 + (D*w)^2 )
    dxdt1 = w * x1;
    F_Coriolis = 2 * M * w_in * dxdt1
    x2_res = F_Coriolis / (D * w20)
    x_max = F1 / (D * w10)
    Q_soi_high = w10 * M / D
    BW_soi_Hz = D / M * 1/(2*pi)
    
%% Temperature    
    thermalexpansion = 2.6e-6;
    Tmin = -55;
    Tmax = 140;
    Troom = 23;
    i = 0;
    for T = Tmin : (Tmax-Tmin)/100 : Tmax
        i = i + 1;
        h = 8*u;
        L2 = 400*u;
        E = 170e9;
        w = (K_soi * L2^3 / h / (E * 192))^(1/3);
        dT = T - Troom;
        L2 = L2 + thermalexpansion * L2 * dT;
        w = w + thermalexpansion * w * dT;
        h = h + thermalexpansion * h * dT;
        K_T = E * 192 * h * w^3 / L2^3;
        V_drive = 15;
        F0 = 2 * 1/2 * N * eps * h / g * V_drive^2;
        %x(i) = F0 / sqrt( (K_T - M*omega^2)^2 + (D*omega)^2 );
        omega_drive(i) = sqrt(K_soi/M);
        omega_T(i) = sqrt(K_T/M);
        TT(i) = T;
    end
    figure(9); clf; grid on; hold on;    
    plot(TT,omega_T/1e3/(2*pi)); xlabel('Temperature [C]'); ylabel('Resonance frequency [kHz]'); title('Comb-IMU.   Input disturbance: Temperature');    
    figure(10); clf; grid on; hold on;    
    plot(TT,(omega_T-omega_drive)/(2*pi)); xlabel('Temperature [C]'); ylabel('Change in resonance frequency [Hz]'); title('Comb-IMU.   Input disturbance: Temperature');    
    
    K_soi
    