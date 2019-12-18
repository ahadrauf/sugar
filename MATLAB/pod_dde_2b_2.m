function maxCurrrent = pod_dde_2b(N_drive,M,D, h,g,V,K, N_sense,N_total)
format compact
%calls: none
%called by: none 
%by: jvclark
%This solves a nonlinear PODMEMS delay differential equation, where
%   Fn = Mechanical and electrical noises
%   Fdr = Excitation
%   FM = Me feedback due to delayed acceleration
%   FD = De feedback due to delayed velocity
%   FK = Ke feedback due to delayed displacement 
%   Fk = ke feedback due to delayed displacement cubed

%Time delays
    tauK = 50e-9; %Time delay of FK
    tauD = tauK; %Delay of FD
    tauM = tauK; %Delay of FM
    tau = [tauK, tauD, tauM];
%Initial conditions
    history = [0, 0]; %Displacement and velocity prior to t=0
    t_initial = 0; %Start time
    t_final   = 1e-3; %Final time
    %t_final   = 7e-3; %Final time
    tspan = [t_initial t_final]; %Simulation time span
%Discontinuity jumps, when the dde coefficient change values
    options = ddeset('jumps',[2,3,4,5,6]*1e-3);
    %options = ddeset('reltol',1e-3,'abstol',1e-5); %AbsTol 1e-6 RelTol 1e-3 
%Solution
    tic %Start timer
    solution = dde23( @PODMEMS, tau, history, tspan, options, tau,N_drive,M,D,h,g,V,K);
    MINUTES = toc/60; %Show how many minutes it took    

%Calculate velocity and current
   count=0; 
   for ii=1:length(solution.y)-1 %number of x after meshed 
       count=count+1;
       %use difference of y when difference is very small to do
       %approximation of the differential.
       Velocity(ii)=(solution.y(1,ii+1)-solution.y(1,ii))/(solution.x(ii+1)-solution.x(ii));
       Current_1(ii)=100*10/2*8.854e-12*Velocity(ii);
       current(count)=Current_1(ii)*(N_sense/N_total);
       time(count)=solution.x(ii);
   end
    maxCurrrent = abs(max(max(Current_1))); %Get the amplitude of the alternating current.
        
%% Plot displacement vs time
%     figure(1); clf; grid on; hold on;
%     title('PODMEMS DDE2, displacement vs time'); xlabel('Time [sec]'); ylabel('Displacement [m]');
%     plot(solution.x, solution.y(1,:),'b'); %Displacement
    
%% Plot current vs time
%     figure(2);grid on; hold on
%     title('PODMEMS DDE23,current vs time'); xlabel('Time [sec]'); ylabel('Current [A]');
%     plot(time,current(1,:),'g');%current
   

%% PODMEMS_1

function dx_dt = PODMEMS(t, x, x_delay, tau, N,M, D, h,g,V,K)
%DDE of PODMEMS
%   t = current time 
%   x = [disp(t); vel(t)] 
%   x_delay = [ disp{tauK), disp{tauD), disp{tauM)
%                vel{tauK),  vel{tauD),  vel{tauM) ]
%   tau = [tauK, tauD, tauM]

%units
    m = 1e-3; %milli
    u = 1e-6; %micro
%Current state
    x1 = x(1);
    x2 = x(2);
%Delayed displacements x1
    x1_tauK = x_delay(1,1);
    x1_tauD = x_delay(1,2);
    x1_tauM = x_delay(1,3);
%Delayed velocities x2
    x2_tauK = x_delay(2,1);
    x2_tauD = x_delay(2,2);
    x2_tauM = x_delay(2,3);
%MEMS structural parameters
%     M = 1.5e-8; %kg
%     D = 6.47e-7; %Ns/m
%     K = 27; %N/m
    k = pi * pi * 170e9 * 2e-6 * 20e-6 / (294.7e-6)^3 / 64; %4e10 [N/m^3]
%Comb params
   % N = 100; %Number of fingers
    e0 = 8.854e-12; %Permittivity of the medium
%     h = 20*u; %Layer thickness
%     g = 20*u; %Gap between comb fingers
%     V = 40; %Voltage 
%Comb force is a square wave
    t0=0; t1=1*m; t2=2*m; t3=3*m; t4=4*m; t5=5*m; t6=6*m; t7=7*m; t8=8*m; T=0.5*m;
    if (t0<=t & t<t0+T)||(t1<=t & t<t1+T)||... %(t2<=t & t<t2+T)||(t3<=t & t<t3+T)||
            (t4<=t & t<t4+T)||(t5<=t & t<t5+T)%||(t6<=t & t<t6+T)
        Fdr = (1/2) * N * 2 * e0 * h * V^2 / g; %Comb drive force
    elseif (t2<=t & t<t2+T)
        Fdr =(1/2) * N * 2 * e0 * h * V^2 / g; %Comb drive force
    elseif (t3<=t & t<t3+T)
        Fdr = (1/2) * N * 2 * e0 * h * V^2 / g; %Comb drive force
    else 
        Fdr = 0; %Off
    end
%Noise
    a = -1; b = 1; %Random interval
    noise_factor = 1e-8; %
    Fn = (a + (b-a).*rand(1,1)) * noise_factor; %Noise always on
%Ke De Me 
    if t0<=t & t<t1
        Ke = 0; De = 0; Me = 0; %Purely mechanical
    elseif t1<=t & t<t2 
        Ke = 0*K; De = 0*D; Me = -0.75*M; %Unstable
    elseif t2<=t & t<t3
        Ke = -0.8*K; De = 200*D; Me = 2*M; %Lower K
    elseif t3<=t & t<t4
        %Ke = -0.5*K; De = 0*D; Me = -0.75*M; %Stable, low K
        Ke = 10*K*(4e-3-t)/1e-3; De = 100*D; Me = 0; %High K and D
    elseif t4<=t & t<t5
        Ke = 0; De = sqrt(4*M*K)/10 - D; Me = 0; %Under-damped
    elseif t5<=t & t<t6
        Ke = 0; De = sqrt(4*M*K) - D; Me = 0; %Over-damped
    elseif t6<=t & t<t7
        Ke = 0; De = 0; Me = -0.93*M; %self-osc
        %elseif t7<=t & t<t8
    else
        Ke = 0; De = 0; Me = 0; %Off
    end
%Feedback forces    
    FK = Ke * x1_tauK; %Feedback force for Ke
    FD = De * x2_tauD; %Feedback force for De
    tauD = tau(2); %Delay 
    tauM = tau(3); %Delay
     x3_tauM = (x2-x2_tauM)/tauM; %Delayed accel at t-tauM
%Net forces + disturbances
    FM = Me * x3_tauM; %Feedback force for Me
    ke = k * 10; %Electrical nonlinear stiffness
    Fk = ke * x1_tauK^3; %Electrical nonlinear stiffness force
    Fnet = 0*Fn + Fdr - FK - Fk - FD - FM; %Totality of all forces
%Delay Differential Equation of PODMEMS
    dx1_dt = x2; %Velocity
    dx2_dt = M \ (Fnet - D*x2 - K*x1 - k*x1^3); %Acceleration
    dx_dt = [dx1_dt; dx2_dt]; %Output
%EOF
 