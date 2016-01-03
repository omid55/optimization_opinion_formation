function fvals = DistributedProjectedSubgrad(loss_funcs, grad_funcs, A, x_init, ...
				     T, stepsize, R)
% DISTRIBUTEDPROJECTEDSUBGRAD Minimizes the given loss function over a graph ...
%using the Distributed Subgradient Projection algorithm from Ram et al, ICASSP, 2009.
%
% fvals = DistributedMinimize(loss_funcs, grad_funcs, A, x_init, T,
%                             stepsize, R)
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
%   The function values of y(t), the projection of the average of the
%   gradients at time t, are returned in fvals.

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

fvals = zeros(T + 1, 1);
for t = 1:T
  % Compute losses for each objective function
  for j = 1:n
    fvals(t) = fvals(t) + loss_funcs{j}(y);
  end
  fvals(t) = fvals(t) / n;
  % Now compute local gradients and update the Z vector
  for j = 1:n
    G(:, j) = grad_funcs{j}(X(:, j));
  end
  Z = X - stepsize*G/sqrt(T);
  parameter_norms = sqrt(sum(Z.^2));
  if (any(parameter_norms > R))
    large_indices = (parameter_norms > R);
    Z(:, large_indices) = R * Z(:, large_indices) ./ ...
	repmat(parameter_norms(large_indices), d, 1);
  end

  X = P*Z;
  y = mean(Z, 2);
  if (norm(y) > R)
      error('norm of y should be smaller than R\n');
  end
  if (mod(t, 100) == 0)
    fprintf(1, 'Distributed Minimization Iteration %d, objval %2.3f\n', ...
            t, fvals(t));
  end
end

for j = 1:n
  fvals(T+1) = fvals(T+1) + loss_funcs{j}(y);
end
fvals(T + 1) = fvals(T + 1) / n;
