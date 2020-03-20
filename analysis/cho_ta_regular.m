% Code to test a version of the simulation algorithm with no scaling.
% For test purposes only.

function [T,Q,M,D,K,l,u,p,il,iu,ip] = cho_ta_regular(net,tspan,varargin) %net,tspan,q0

% -- Initialize
numvararg=length(varargin); %number of varargs
T=0; 

% -- Setup system
global M D K;

M = assemble_system(net, 'M');
D = assemble_system(net, 'D');
K = assemble_system(net, 'K');
k = length(K);
[l,u,p]=lu(M);
ip=pinv(p);
il=pinv(l);
iu=pinv(u);
invMD=pinv(D)*M;%M\D;
invMK=pinv(K)*M;%M\K;
%disp(M);
%disp(D);
%disp(K);
%return;

% -- Set initial conditions

if numvararg >=1 
  q0 = varargin{1};
  if isempty(q0) %if q0 is [].
    q0 = zeros(2*k,1); %set it to 0.
  end
else
  q0 = zeros(2*k,1); %set it to 0.
end

% -- Integrate

%options = odeset('Mass', 'M');
%[T,Q] = ode23s('ta_ode', tspan, q0, options, net, C, G, d1, d2);
[T,Q] = ode45('ta_ode_regular', tspan, q0, [], net, M, D, K,l,u,p,il,iu,ip,invMD,invMK);
