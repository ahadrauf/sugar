% Find SS response of linearized system about the equilibrium point q0

% Note: Besides needing to document and test this one, I *seriously*
%       need to sanity check the influence and output matrices.

function find_ss(net, q0, in_node, in_var, out_node, out_var)


% -- Determine influence matrix

B = zeros(2*net.dof,1);
in_var_id = lookup_coord(net, in_node, in_var)
if (in_var_id ~= 0)
  B(in_var_id) = 1;
end

% -- Determine output matrix

C = zeros(1,2*net.dof);
out_var_id = lookup_coord(net, out_node, out_var)
if (out_var_id ~= 0)
  C(out_var_id) = 1;
end


% -- Assemble scaled linear system

[C1,G1, D1,d2] = setup_linear(net,q0);

% -- Rescale the influence and output matrices

B = B./d2;
B = D1*B;
C = C.*(d2');


% -- Put together state space model (doing straight tf would actually
%    probably be more stable in this case...)

A = -C1\G1;
B = C1\B;
sys = ss(A,B,C,0);


% -- Draw Bode plot

bode(sys);
