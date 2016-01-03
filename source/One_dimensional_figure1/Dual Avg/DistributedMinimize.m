function [centralized_fvals local_fvals x_avg x_avg_time] = DistributedMinimize(loss_funcs, grad_funcs, A, x_init, T, stepsize, R, steptype)
% DISTRIBUTEDMINIMIZE Minimizes the given loss function over a graph.
%
% [centralized_fvals local_fvals x_avg] = ...
%     DistributedMinimize(loss_funcs, grad_funcs, A, x_init, T, ...
%                         stepsize, R, steptype)
%   minimizes the sum of the loss functions stored in the cell array
%   loss_funcs. The function assumes that loss_funcs and grad_funcs are
%   cell arrays of function pointers, so we minimize
%
%    sum_i (loss_funcs{i}(x))
%
%   subject to an L2 constraint on the parameter vectors, that is, norm(x) <=
%   R. The method used is the standard distributed lazy projection algorithm
%   started from initial starting point 'x_init', run for T iterations using
%   the given stepsize, and adjacency matrix A for the graph. A is assumed
%   to include self-loops.
%
%   The function values of y(t), the projection of the average of the gradients
%   at time t, are returned in centralized_fvals, and the function values of
%   x_1(t), the parameters at node 1, are returned in local_fvals. x_avg is
%   the average of all x values at node 1.
%
%   If steptype == 1, uses fixed step rate. If steptype == 2, uses 1/sqrt(t)
%   stepping.

n = size(A, 1);
if (n ~= numel(loss_funcs) || n ~= numel(grad_funcs))
  error('Wrong number of loss or gradient functions');
end
% First, construct the transition matrix.
P = ConstructTransitionFromAdjacency(A);

d = size(x_init, 1);

% Parameter vectors, local gradient vectors, average gradient collections
X = repmat(x_init, 1, n);
G = zeros(d, n);
Z = zeros(d, n);

y = x_init;

x_avg = zeros(d, 1);
x_avg_time = zeros(d, T);

centralized_fvals = zeros(T + 1, 1);
local_fvals = zeros(T + 1, 1);
for t = 1:T
  % Compute losses for each objective function, for both centralized and local.
  for j = 1:n
    centralized_fvals(t) = centralized_fvals(t) + loss_funcs{j}(y);
    local_fvals(t) = local_fvals(t) + loss_funcs{j}(X(:, 1));
  end
  centralized_fvals(t) = centralized_fvals(t) / n;
  local_fvals(t) = local_fvals(t) / n;
  % Now compute local gradients and update the Z vector
  for j = 1:n
    G(:, j) = grad_funcs{j}(X(:, j));
  end
  Z = Z*P + G;
  % Project to 2-norm ball
  X = -stepsize * Z / sqrt(t);
  parameter_norms = sqrt(sum(X.^2));
  if (any(parameter_norms > R))
    large_indices = (parameter_norms > R);
    X(:, large_indices) = R * X(:, large_indices) ./ ...
	repmat(parameter_norms(large_indices), d, 1);
  end
  y = -stepsize * mean(Z, 2) / sqrt(t);
  if (norm(y) > R)
    y = R * y / norm(y);
  end
  if (mod(t, 100) == 0)
    fprintf(1, ['Distributed Minimization Iteration %d, centralized' ...
		'objval %2.3f, local objval %2.3f\n'], ...
            t, centralized_fvals(t), local_fvals(t));
  end
  x_avg = x_avg + X(:, 1);
  x_avg_time(:, t) = x_avg / t;
end

for j = 1:n
  centralized_fvals(T+1) = centralized_fvals(T+1) + loss_funcs{j}(y);
  local_fvals(T + 1) = local_fvals(T + 1) + loss_funcs{j}(X(:, 1));
end
centralized_fvals(T + 1) = centralized_fvals(T + 1) / n;
local_fvals(T + 1) = local_fvals(T + 1) / n;
x_avg = x_avg / T;
