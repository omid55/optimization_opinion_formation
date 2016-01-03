function [opt_val B obj_funcs grad_funcs] = MakeSVMData(n, d, noise, G, R)
% MAKESVMDATA Creates objective functions to minimize in distributed fashion
%
% [opt_val B obj_funcs grad_funcs] = MakeSVMData(n, d, noise, G, R) constructs
%   n different functions, each of which is a hinge loss. The functions take
%   d-dimensional vectors as arguments, and G is the the Lipschitz constant of
%   the functions (i.e. norm(grad_funcs{i}(x)) <= G). If G is unspecified, 1 is
%   used. obj_funcs and grad_funcs are cell arrays of size n, each entry of
%   which contains a function pointer to either a gradient or objective
%   function.
%
%   Specifically, the data is constructed as an n-by-d random matrix B whose
%   rows are all norm G. A vector w is chosen randomly, and labels y =
%   sign(B*w) are set. The noise parameter (.05 if unspecified) controls the
%   number of elements of y that are flipped. A randomly selected set of noise
%   * numel(y) elements of y have their signs flipped.
%
%   Each objective function is then
%
%   max(1 - y(i) * B(i, :) * x, 0)
%
%   and the global objective is their average. The opt_val returned is the
%   optimal value of the minimization problem
%
%   min.  sum_i max(1 - y(i) * B(i, :) * x, 0)
%   s.t.  norm(x) <= R
%
%   If R is unspecified, uses R = 5.
%
% Author: John Duchi (jduchi@cs.berkeley.edu)

if (nargin <= 2)
  noise = .1;
end
if (nargin <= 3)
  G = 1;
end
if (nargin <= 4)
  R = 5;
end

B = randn(n, d);
B = G * B ./ sqrt(repmat(sum(B.^2, 2), 1, d));
% We use the noise level to make separation more difficult
y = zeros(n, 1);
x = randn(d, 1);
if (noise < 1)
  for ii = 1:n
    if (rand > noise)
      y(ii) = sign(B(ii, :) * x);
    else
      y(ii) = sign(1 - 2 * rand);
    end
  end
else
  y = ones(n, 1);
end
% y = sign(B * randn(d, 1));
% flip_inds = rand(n, 1) < noise;
% y(flip_inds) = -y(flip_inds);
B = repmat(y, 1, d) .* B;
fprintf(1, 'Solving problem with cvx... ');
cvx_quiet(true);
cvx_begin
  variable x(d);
  minimize(sum(max(1 - B*x, 0)));
  norm(x) <= R;
cvx_end
fprintf(1, 'done.\n');

opt_val = sum(max(1 - B * x, 0)) / n;

if (nargout > 2)
  obj_funcs = cell(n, 1);
  grad_funcs = cell(n, 1);
  for ii = 1:n
    obj_funcs{ii} = @(x) HingeLoss(B(ii, :), 1, x);
    grad_funcs{ii} = @(x) HingeGrad(B(ii, :), 1, x);
  end
end
