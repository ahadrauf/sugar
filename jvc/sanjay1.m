clear all;

% Solves for delta-V where V(x,y) = V1*x/d + delta-V
% Assumes periodic BC; i.e. V(-L/2,y) = V(L/2,y)

L = 20;      % Length of gap
d = 1;       % Average gap value
dd = 0.01*d; % Perturbation of gap must be small for validity of result
m = 10;      % wave number for perturbation
V1 = 1;      % Upper voltage ; lower voltage = 0;
terms = 6;   % Number of Fourier terms

% Build RHS
rhs = zeros(terms,1);
rhs(1) = V1*dd/d;

% Build Matrix
for jj = 1:terms
    A(jj,jj) = sinh(2*pi*jj*m*d/L);
    if jj+1 <= terms
        A(jj,jj+1) = pi*m*(jj+1)*dd*cosh(2*pi*m*(jj+1)*d/L)/L;
    end
    if jj-1 >= 1
        A(jj,jj-1) = pi*m*(jj-1)*dd*cosh(2*pi*m*(jj-1)*d/L)/L;
    end
end

% Solve for coefficients
C = A\rhs;

% Plot
x = linspace(-L/2,L/2,100);
y = linspace(0,d+dd,100);

V = zeros(length(x),length(y));
dV = zeros(length(x),length(y));

for ii=1:length(x)
    for kk=1:length(y)
        V(ii,kk)  = V1*y(kk)/d;
        dV(ii,kk) = 0;
        for jj = 1:terms
            dVx       = C(jj)*cos(2*pi*jj*m*x(ii)/L)*sinh(2*pi*jj*m*y(kk)/L);
            dV(ii,kk) = dV(ii,kk) + dVx;
            V(ii,kk)  = V(ii,kk) + dVx;
        end
    end
end


% Warning: Plots some points outside the domain!
surf(x,y,V')
title('Total potential');

figure
surf(x,y,dV')
title('Potential variation due to roughness')
        
