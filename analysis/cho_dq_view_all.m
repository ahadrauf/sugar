% function [] = cho_dq_view_all(q, net)
% This function provides a comprehensive view into the state output of all cho_* functions
% Examples:
% 	net = cho_load("cantilever.net");
% 	[q] = cho_dc(net);
% 	cho_dq_view_all(q, net);
%
%	[t, q] = cho_ta(net, [0, 1e-4]);
% 	cho_dq_view_all(q(end,:), net)
%
%Example output (for gap_closing_actuator.net):
%1 = A, e, 80
%2 = b, e, 80
%3 = b, x, 0
%4 = b, y, -4.65573e-12
%5 = b, rz, -8.939e-07
%6 = c, e, 80
%7 = c, x, 0
%8 = c, y, -1.63044e-10
%9 = c, rz, -2.71895e-06
%10 = D, e, 0
%11 = D, x, 0
%12 = D, y, 1.63044e-10
%13 = D, rz, -2.71895e-06
%14 = E, e, 0
%15 = E, x, 0
%16 = E, y, 4.65573e-12
%17 = E, rz, -8.939e-07
%18 = G, i, 0



function [] = cho_dq_view_all(q, net)
    for i=1:length(net.vars)
        if ~strcmp('g', net.vars(i).type)
            owner = net.vars(i).owner;
            if owner <= length(net.nodes)
                node_name = net.nodes(owner).name;
            else
                node_name = "N/A";
            end
            var_name = net.vars(i).name;
            var_value = q(i);
            fprintf("%d = %s, %s, %g\n", i, node_name, var_name, var_value);
        end
    end
end