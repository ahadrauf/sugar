echo off
%SUGAR Tutorial. Electrical-mechanical analysis of beam gap structure.

%   Created by Ningning Zhou on 02/13/2001;
%   Updated: dbindel, 9/8/2001

echo on; clc
%                     SUGAR Tutorial - Static Analysis
%##########################################################################
%   Electrical-mechanical analysis of beam gap structure.%
%##########################################################################

%pause % Strike any key to display netlist "beamgap2e.net"; 
%
%type beamgap2e.net
%
%pause % Strike any key to display process file "mumps.net"; 
%
%type mumps.net
%
%pause % Strike any key to load netlist and display the original structure;5
%
net = cho_load('cantilever.net');
%
%figure;
%cho_display(net);

%pause % Strike any key to continue
%clc

%########################################################################
%
%   Static Analysis;
%
%#######################################################################

%pause % Strike any key to analyze and display deflection at 10V

% [Q] = cho_dc(net);
V = 40;
p.V = V;
net = cho_load('beamgap2e.net');
[T,Q] = cho_ta_regular(net,[0 1e-4]);   % Simulate 1 ms behavior 
dy = cho_dq_view(Q, net, 'c', 'y'); % Get the y component at c 
plot(T, dy);    % plot how it moves over t

%figure;
%cho_display(net,dq);
% title('Deflected structure at V=10v');

%pause % Strike any key to exit
echo off
disp('End of demo')

% net = cho_load('cantilever.net');
% [T,Q] = cho_ta(net,[0 1e-3]);   % Simulate 1 ms behavior 
% dy = cho_dq_view(Q, net, 'tip', 'y'); % Get the y component at c 
% plot(T, dy);    % plot how it moves over t
