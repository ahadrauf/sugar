% Simulate the nonlinear system over some time period.
%
% Inputs:
%   net   - netlist structure
%   tspan - two-element vector [tstart tend] indicating the start
%           and end times for the simulation.
%   q0    - (Optional) initial state at time = tstart.
%           If q0 is not provided, the default is zero.
%
% Output:
%   T     - time points where the solution was sampled
%   Q     - array of state vectors sampled at the times in T
%           (i.e. Q(i,:) is the state vector at time T(i))

function [T,Q,C,G] = cho_ta(net,tspan,q0)


% -- Setup system

[C,G,d1,d2] = setup_linear(net, [], 1);

% -- Set initial conditions if unspecified

if nargin < 3
  q0 = zeros(length(d1),1); 
end

% -- Integrate

options = odeset('Mass', 'M', 'Jacobian', 'on');

[T,Q] = ode23s('ta_ode', tspan, q0, options, net, C, G, d1, d2);
Q = Q*diag(d2);
