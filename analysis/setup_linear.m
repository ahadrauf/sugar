% Setup scaled linearized system about an equilibrium point.
% Note that the system is assumed to be autonomous.
%
% Note also that 
%  Chat^-1 Ghat = (D2^-1 C^-1 D1^-1) (D1 G^-1 D2)
%               = D2^-1 (C^-1 G) D2
% ie Chat\Ghat is similar to C\G.  Consequently, any characterization
% involving the eigenvalues of 
%  Chat zhat' + Ghat zhat = 0 
% will carry over to 
%  C z' + G z = 0
% though the eigenvectors of the systems will differ by a diagonal
% scaling.
%
% Inputs:
%   net   - netlist structure 
%   q0    - equilibrium state
%   no_F  - (optional) If true, form the first-order form of the
%           linear system given by setting F = 0.  In this case,
%           q0 is ignored.
%
% Output:
%   Chat, Ghat - 
%            matrices in the first-order linearization
%             Chat' [zhat; zhat']' + Ghat [zhat; zhat'] = 0
%            where z = D2*zhat = q-q0
%   d1, d2 - diagonals of diagonal scaling matrices D1 and D2
%            used to equilibrate the system (Chat = D1*C*D2,
%            Ghat = D1*G*D2)

function [C,G,D1,d2] = setup_linear(net,q0,no_F);


% -- Set default argument

if (nargin <= 2)
  no_F = 0;
end


% -- Setup system

M = assemble_system(net, 'M');
D = assemble_system(net, 'D');
K = assemble_system(net, 'K');

k = length(K);
if (no_F == 0)
   K = K-assemble_dFdx(net, q0, 0);
   d = -assemble_dFdxdot(net, q0, 0);
else
   d = zeros(k);
end

% -- Convert to first order

W = eye(k);
o = zeros(k);

C = [D M; W o];
G = [K d; o -W];

% -- Compute scalings

for i=1:(2*k)
  nC = norm(C(i,:));
  if nC == 0
    d1(i) = 1;
  else
    d1(i) = 1/nC;
  end
end
d1 = d1';

D1 = diag(d1);
C = D1*C;
G = D1*G; 

for i=1:(2*k)
  %d2(i) = 1/norm(G(:,i));
  d2(i) = 1;
end
d2 = d2';

D2 = diag(d2);
C = C*D2;
G = G*D2;


