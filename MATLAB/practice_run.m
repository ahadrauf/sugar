% '////////////////////////////////'
clear; clc;
tic % start clock

V = 40;
%q0 = [-10.0000  -10.0000         0   -0.0000   -0.0033  -10.0000         0   -0.0000   -0.0033         0]';
%q0 = [20.0000   20.0000         0   -0.0000   -0.0000         0   -0.0000   -0.0000         0    0.0000   -0.0000         0    0.0000   -0.0000         0 0]';
%q0 = [20.0000   20.0000         0    0.0000    0.0014         0    0.0000    0.0062   20.0000         0         0         0         0         0         0]';
%q0 = [20.0000   20.0000         0   -0.0000   -0.0000   20.0000         0   -0.0000   -0.0000         0         0    0.0000   -0.0000         0         0 0 0 0]';
q0 = [4.0000e+01   4.0000e+01            0  -5.6129e-12  -1.0809e-06   4.0000e+01            0  -2.2803e-10  -3.5176e-06            0            0   2.2803e-10 -3.5176e-06            0            0   5.6129e-12  -1.0809e-06            0]';
q0 = [q0; zeros(length(q0), 1)];

p.V = V;
% net = cho_load('beamgap2e.net', p);
net = cho_load('practice_structure.net', p);
%net = cho_load('cantilever.net');
[q] = cho_dc(net);
%[T,Q] = cho_ta(net,[0 1e-5], q0);   % Simulate 1 ms behavior
%dy = cho_dq_view(Q, net, 'c', 'y'); % Get the y component at c 
%plot(T, dy);    % plot how it moves over t

%figure(1);
%cho_display(net);

%[q] = cho_dc(net);


%figure();
%cho_display(net, q);



%%%%%% prototype
q0 = [cho_dq_view(q, net, 'c', 'x'), cho_dq_view(q, net, 'c', 'y'), 0, 0, 0, cho_dq_view(q, net, 'c', 'rz')]';
q0 = [q0; zeros(length(q0), 1)];
net = cho_load('cantilever.net');
%[q] = cho_dc(net);
[T,Q] = cho_ta(net,[0 0.25e-5], q0);   % Simulate 1 ms behavior
dy = cho_dq_view(Q, net, 'tip', 'y'); % Get the y component at c 
plot(T, dy);    % plot how it moves over t

return

p.V = 0;
net = cho_load('beamgap2e.net', p);
[q] = cho_dc(net, q0);
figure(3);
cho_display(net, q);



toc % stop clock
