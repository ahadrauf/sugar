% '////////////////////////////////'
clear; clc;
tic % start clock

voltages = linspace(40, 80, 5);
for i=1:length(voltages)
    [temp1, temp2] = ta_gca_helper(voltages(i));
    T{i} = temp1;
    dy{i} = temp2;
end

figure(2);
hold on;
for i=1:length(voltages)
    plot(T{i}, dy{i}, 'DisplayName', sprintf('V = %d', voltages(i)));    % plot how it moves over t
end
yline(max(dy{1}), 'DisplayName', 'Gap Stop');
legend

toc % stop clock

function [T, dy] = ta_gca_helper(V)
%% Calculate equilibrium point of pull-in
    p.V = V;
    net = cho_load('gap_closing_actuator.net', p);
    [q] = cho_dc(net);

    %figure(1);
    %cho_display(net); % Figure before displacement
    %cho_display(net, q); % Figure after displacement
%% Turn equilibrium point of pull-in to initial position of pull-out
    q0 = [cho_dq_view(q, net, 'c', 'x'), cho_dq_view(q, net, 'c', 'y'), 0, 0, 0, cho_dq_view(q, net, 'c', 'rz')]';
%     q0 = [cho_dq_view(q, net, 'c', 'x'), -1.0282e-10, 0, 0, 0, cho_dq_view(q, net, 'c', 'rz')]';
    q0 = [q0; zeros(length(q0), 1)];
    net = cho_load('gap_closing_actuator_mechanical_analog.net');
    [T,Q] = cho_ta(net,[0 0.5e-6], q0);   % Simulate 2.5 us behavior (pull-out)
    dy = cho_dq_view(Q, net, 'tip', 'y'); % Get the y component at the gap closing actuator finger tip
end