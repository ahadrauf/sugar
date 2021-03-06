% This function runs a DC simulation over an angled arm, implemented in a
% 2-mask SOI process. It assesses the displacement gain  (delta X to delta
% Y) and the force gain (delta X to F_x and F_y to F_x) on the shuttle.
% It's divided into three portions: the first runs a straightforward sweep
% of F_y to F_x with fixed L_arm and alpha. The second part fixes F_y
% and L_arm and varies alpha. The third part fixes F_x and alpha and
% varies L_arm.
%
% This simulation loads the subnet_angled_arm.net subnet, located in the
% demo/ folder.
% 
% As is, this simulation takes ~10-15 minutes to run. This time can be
% reduced by reducing the number of L_arm and alpha values that are
% tested. These can be adjusted with the num_vals variables (or by directly
% manipulating alpha_vals and L_vals).

%force output over 2um displacement
clear;clc;

% define constants and variables
E=169e9;
alpha=65;
Larm=125.98e-6;
warm=3e-6;
% t=40e-6;
t = 550e-6;
Iarm=warm^3*t/12;
Kphi=E*Iarm/Larm;
Wb=0.1e-6;
% Fy=(100:300:2000)*1e-6;
V = [30, 40, 50, 60, 80, 100];
Fx=0;

% define parameter structure for netlist
param.Fx=Fx;
% param.Fy=Fy(end);
param.V = V(end);
param.Wb=Wb;
param.Larm=Larm;
param.warm=warm;
param.alpha=alpha*pi/180;
param.V = 60;
% take first step with netlist
net = cho_load('subnet_angled_arm.m',param);
dq1 = cho_dc(net);

% Display figure
figure(1); cho_display(net);
figure(2); cho_display(net,dq1);
return;

% Run the first sweep F_x vs. F_y with fixed L_arm and alpha
dy1 = dqval(net,dq1,'n2','y');
dy2 = dqval(net,dq1,'n1','y');
dx1 = dqval(net,dq1,'p1','x');
dx2 = dqval(net,dq1,'nmid','x');
dstop=0e-6;
% dstop = -3.833e-6 * tand(65);
for i=1:length(V)
    Fx=0;
    param.Fx=Fx;
%     param.Fy=Fy(i);
    param.V = V(i);
    dq = dq1;
    while dx2>dstop
        net = cho_load('subnet_angled_arm.m',param);
        dq = cho_dc(net);
        dy1 = dqval(net,dq,'n1','y');
        dy2 = dqval(net,dq,'n2','y');
        dx1 = dqval(net,dq,'p1','x');
        dx2 = dqval(net,dq,'nmid','x');
        
        if Fx == 0
            dx2_orig = dx2;
        end
        
        Fx=Fx+10e-6;
        param.Fx=Fx;
    end
    figure; cho_display(net,dq);
    Fout(i)=Fx;
%     Foutcalc(i) = Fy/tand(alpha) - 3*E*Iarm/Larm^3 * dx2_orig/sind(alpha)^2;
    dx1=1;
    dx2=1;
    dy1=-1;
    dy2=1;
end
% Foutcalc=2*(Fy/tand(alpha));

figure(3);
hold on;
% plot(Fout*1e6,Fy*1e6, 'LineWidth', 4);
plot(V, Fout*1e6, 'LineWidth', 4);
% plot(Foutcalc*1e6,Fy*1e6, 'LineWidth', 4);
xlabel('V'); ylabel('F_y (uN)');
title('V vs. F_x');
% legend('SUGAR', 'Analytical');

return
% Run the second sweep, alpha vs F_x, with fixed F_y and L_arm
num_vals = 5;
alpha_vals = linspace(60, 80, num_vals);
for L=1:length(alpha_vals)
    alpha = alpha_vals(L);
    param.alpha = alpha_vals(L) * pi/180;
    
    Fy=1000e-6;
    Fx=0;
    param.Fy=Fy;
    param.Fx=Fx;
    net = cho_load('subnet_angled_arm.net',param);
    dq = cho_dc(net);
    dx = dqval(net,dq,'p1','x');
    i=1;
    while dx>5e-6
        param.Fx=param.Fx+1e-6;
        net = cho_load('subnet_angled_arm.net',param);
        dq = cho_dc(net);
        dx = dqval(net,dq,'p1','x');
        dy = dqval(net,dq,'n2','y');
    end
    dx1A(L,1)=dx;
    dy1A(L,1)=dy;
    while dx1A(L,i)>0
        i=i+1;
        dx1A(L,i) = dqval(net,dq,'p1','x');
        FxoutA(L,i)=param.Fx;
        param.Fx=param.Fx+1e-6;
        net = cho_load('subnet_angled_arm.net',param);
        dq = cho_dc(net);
        dy1A(L,i) = dqval(net,dq,'n2','y');
    end
    FxcalcA(L,:)=2*(Fy/tand(alpha)-3*Kphi*dx1A(L,:)/(Larm^2*sind(alpha)^2));
    dxcalcA(L,:)=dy1A(L,:)*tand(alpha);
end

figure(4);
hold on;
temp1 = dx1A;
temp2 = dy1A;
temp1(temp1 == 0 | temp2 == 0) = NaN;
temp2(temp1 == 0 | temp2 == 0) = NaN;
for L=1:length(alpha_vals)
    plot(temp1(L,:)*1e6,temp2(L,:)*1e6, 'LineWidth', 4, 'DisplayName',sprintf('alpha = %g�', alpha_vals(L)));
    %plot(dxcalcA(L,:)*1e6,dy1A(L,:)*1e6, 'LineWidth', 4);  % Analytical
end
legend
xlabel('\Deltax (um)'); ylabel('\Deltay (um)');
title('\Deltax vs. \Deltay');

figure(5);
hold on;
temp1 = dx1A;
temp2 = FxoutA; 
temp1(temp1 == 0 | temp2 == 0) = NaN;
temp2(temp1 == 0 | temp2 == 0) = NaN;
for L=1:length(alpha_vals)
    plot(temp1(L,:)*1e6,temp2(L,:)*1e6, 'LineWidth', 4, 'DisplayName',sprintf('alpha = %g�', alpha_vals(L)));
    %plot(dx1A(L,:)*1e6,FxcalcA(L,:)*1e6, 'LineWidth', 4);  % Analytical
end
legend;
xlabel('\Deltax (um)'); ylabel('F_x (uN)');
title('\Deltax vs. F_x');

% Runs the third sweep, L_arm vs. F_x, with fixed F_y and alpha
alpha = 65;
param.alpha = 65 * pi/180;
num_vals = 5;
Larm_vals = linspace(180e-6, 1000e-6, num_vals);
for L=1:length(Larm_vals)
    Larm = Larm_vals(L);
    param.Larm = Larm_vals(L); % vary L_arm
    
    Fy=100e-6;
    Fx=0;
    param.Fy=Fy;
    param.Fx=Fx;
    net = cho_load('subnet_angled_arm.net',param);
    dq = cho_dc(net);
    dx = dqval(net,dq,'p1','x');
    i=1;
    while dx>5e-6
        param.Fx=param.Fx+1e-6;
        net = cho_load('subnet_angled_arm.net',param);
        dq = cho_dc(net);
        dx = dqval(net,dq,'p1','x');
        dy = dqval(net,dq,'n2','y');
    end
    dx1L(L,1)=dx;
    dy1L(L,1)=dy;
    while dx1L(L,i)>0
        i=i+1;
        dx1L(L,i) = dqval(net,dq,'p1','x');
        FxoutL(L,i)=param.Fx;
        param.Fx=param.Fx+1e-6;
        net = cho_load('subnet_angled_arm.net',param);
        dq = cho_dc(net);
        dy1L(L,i) = dqval(net,dq,'n2','y');
%         i=i+1;
    end
    FxcalcL(L,:)=2*(Fy/tand(alpha)-3*Kphi*dx1L(L,:)/(Larm^2*sind(alpha)^2));
    dxcalcL(L,:)=dy1L(L,:)*tand(alpha);
end

figure(6);
hold on;
temp1 = dx1L;
temp2 = dy1L; 
temp1(temp1 == 0 | temp2 == 0) = NaN;
temp2(temp1 == 0 | temp2 == 0) = NaN;
for L=1:length(Larm_vals)
    plot(temp1(L,:)*1e6,temp2(L,:)*1e6, 'LineWidth', 4, 'DisplayName',sprintf('L = %g um', Larm_vals(L) * 1e6));
    %plot(dxcalcL(L,:)*1e6,dy1L(L,:)*1e6, 'LineWidth', 4);  % Analytical
end
legend
xlabel('\Deltax (um)'); ylabel('\Deltay (um)');
title('\Deltax vs. \Deltay');

figure(7);
hold on;
temp1 = dx1L;
temp2 = FxoutL; 
temp1(temp1 == 0 | temp2 == 0) = NaN;
temp2(temp1 == 0 | temp2 == 0) = NaN;
for L=1:length(Larm_vals)
    plot(temp1(L,:)*1e6,temp2(L,:)*1e6, 'LineWidth', 4, 'DisplayName',sprintf('L = %g um', Larm_vals(L) * 1e6));
    %plot(dx1L(L,:)*1e6,FxcalcL(L,:)*1e6, 'LineWidth', 4);  % Analytical
end
legend;
xlabel('\Deltax (um)'); ylabel('F_x (uN)');
title('\Deltax vs. F_x');
