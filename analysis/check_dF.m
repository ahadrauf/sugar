% Simple function to visually compare the computed Jacobian
% for an element to a difference approximation of the Jacobian.
%
% Inputs:
%   net - name of a netlist structure
%   elt_num - a single element from the element array of a netlist structure
%             (e.g. net.element(1))
%   q - displacement of all coordinates
%   t - time to evaluate F and dF
%   k - size of blocks to compare (for a gap with 6 dof per node, for
%       example, 6 would be a reasonable choice).  Should divide the number
%       of rows/columns in the Jacobian.

function check_dF(net, elt_num, q, t, k)

element = net.elements(elt_num);
model = element.model;
param = element.parameter;
R = element.R;

j = find(element.var_ids ~= 0);
jdx = element.var_ids(j);

x = zeros(2*length(element.var_ids),1);
x(j) = q(jdx);

for j = 1:length(x)
  h = sqrt(eps);
  ej = zeros(length(x),1); 
  ej(j) = 1;
  dFdiff(:,j) = ( feval(model, 'F', R, param, x + h*ej, t) - ...
                  feval(model, 'F', R, param, x - h*ej, t) ) / (2*h);
end

dF = feval(model, 'dFdx', R, param, x, t);

format short e
[r,c] = size(dF)
for i = 1:r/k
  for j = 1:c/k
    is = (i-1)*k + (1:k);
    js = (j-1)*k + (1:k);
    fprintf('Compare %d:%d by %d:%d blocks\n', is(1), is(k), js(1), js(k));
    disp('dF model:');
    disp(dF(is,js))
    disp('dF divided differences:');
    disp(dFdiff(is,js))
    pause
  end
end
