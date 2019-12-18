'////////////////////////////////'
clear all
clear y1 y2 vv v1 v2 ii vv i
V1 = 0;
V2 = 20; 
V3 = 20;
V4 = 0.001; 
V5 = 0.001;
V6 = 0;

tic
i = 0;

%init
    p.V = V1; 
    net = cho_load('net_ucb_contact_zip_1.m',p);
    [q] = cho_dc_q(net);

%for V = [V1:(V2-V1)/100:V2]
for V = [ V1 : (V2-V1)/100 : V2,  V3 : (V4-V3)/100 : V4, V5 : (V6-V5)/1000 : V6, 0, 0]
    i = i + 1;
    p.V = V; 
    net = cho_load('net_ucb_contact_zip_1.m',p);
    [q] = cho_dc_q(net,q);
    %figure(1); cho_display(net,q);pause
    y1(i) = q(lookup_coord(net,'a(1)','y'));
    y2(i) = q(lookup_coord(net,'a(4)','y'));
    v1(i) = q(lookup_coord(net,'a(1)','e'));
    v2(i) = q(lookup_coord(net,'a(4)','e'));
    vv(i) = V;
    ii(i) = i;
end

figure(1);
cho_display(net,q);

figure(2); clf; hold on; grid on;
plot(ii,y1,'b');
plot(ii,y2,'r');
y1=y1(i)
y2=y2(i) 
   
figure(3); clf; hold on; grid on;
plot(ii,v1,'b');
plot(ii,v2,'r');

figure(2);
title('Deflection [node a1 (blue), a2 (red)] vs anchor voltage (node e)');
xlabel('Input anchor voltage of node e [V]');
ylabel('Output deflection of node a1 (blue) & a2 (red) [m]');

figure(3);
title('Gap voltage [node a1 (blue), a2 (red)] vs anchor voltage (node e)');
xlabel('Input anghor voltage of node e [V]');
ylabel('Output gap voltage of node a1 (blue) & a2 (red) [m]');
toc
