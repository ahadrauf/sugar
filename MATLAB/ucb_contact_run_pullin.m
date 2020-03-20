% '////////////////////////////////'
clear y1 y2 vv v1 v2
tic
V1 = 8;
V2 = 8; 
i = 0;

%init
p.V = V1; 
net = cho_load('net_ucb_contact_1b.m',p);
% cho_display(net);
%[T,Q] = cho_ta(net,{0, 1e-3}); % Simulate 1 ms behavior
%dy = cho_dq_view(Q, net, 'a', 'y'); % Get the y component at c, and
%plot(T, dy); % plot how it moves over time

%%%%
[q] = cho_dc_q(net);
cho_display(net, q);
return;


    
for V = V1 : (V2-V1)/1000 : V2
    i = i + 1;
    p.V = V; 
    net = cho_load('net_ucb_contact_1b.m',p);
    [q] = cho_dc_q(net,q);
    y1(i) = q(lookup_coord(net,'a','y'));
    y2(i) = q(lookup_coord(net,'b','y'));
    v1(i) = q(lookup_coord(net,'a','e'));
    v2(i) = q(lookup_coord(net,'b','e'));
    vv(i) = V;
end

figure(1);
cho_display(net,q);

figure(2); clf; hold on; grid on;
plot(vv,y1,'b.-');
plot(vv,y2,'r.-');
Y1 = y1(i)
Y2 = y2(i) 
V1 = v1(i)
V2 = v2(i) 
   
figure(3); clf; hold on; grid on;
plot(vv,v1,'b.-');
plot(vv,v2,'r.-');

figure(2);
title('Deflection [node a1 (blue), a2 (red)] vs anchor voltage (node e)');
xlabel('Input anchor voltage of node e [V]');
ylabel('Output deflection of node a1 (blue) & a2 (red) [m]');

figure(3);
title('Gap voltage [node a1 (blue), a2 (red)] vs anchor voltage (node e)');
xlabel('Input anghor voltage of node e [V]');
ylabel('Output gap voltage of node a1 (blue) & a2 (red) [m]');
toc
