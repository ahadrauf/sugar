% '////////////////////////////////'
clear; clc;
tic % start clock

V = 40;
q0 = [-10.0000  -10.0000         0   -0.0000   -0.0033  -10.0000         0   -0.0000   -0.0033         0]';

p.V = V;
% net = cho_load('beamgap2e.net', p);
net = cho_load('practice_structure.net', p);

%figure(1);
%cho_display(net);

[q] = cho_dc(net);

%figure();
cho_display(net, q);

return

p.V = 0;
net = cho_load('beamgap2e.net', p);
[q] = cho_dc(net, q0);
figure(3);
cho_display(net, q);



toc % stop clock
