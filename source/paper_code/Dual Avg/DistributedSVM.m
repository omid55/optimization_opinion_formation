function [fvals meanFits] = DistributedSVM(B, A, X, T, stepsize, R,optimal)
% DISTRIBUTEDSVM Simulates distributed minimization of hinge loss
%
% fvals = DistributedSVM(B, A, x_init, T, stepsize, R) simulates distributed
%   minimization of the function
%
%   sum(max(1 - B*x, 0))
%
%   over n computers, each of which has one row of B as its objective (i.e. one
%   term in the sum). In the above, B is an n-by-d matrix and x is a
%   d-vector. Thus x_init should be a d-vector. Minimization is performed
%   subject to an L2 constraint on the parameter vector, that is, norm(x) <=
%   R. The method used is standard distributed lazy projection started from
%   x_init. The method is run for T iterations with 1/sqrt(t) decreasing
%   stepsizes, with initial stepsize provide. A is the adjacency matrix of
%   the underlying graph and is assumed to include self-loops.
%
%   The function values of y(t), the projection of the average of the
%   gradients at time t, are returned in fvals.
%
% Author: John Duchi (jduchi@cs.berkeley.edu)
% Reference: Agarwal, Duchi, and Wainwright 2010.

X = X';
n = size(A, 1);
if (n ~= size(B, 1))
  error('Wrong number of loss functions versus network size\n');
end
d = size(X, 1);
if (d ~= size(B, 2))
  error('Wrong initial vector size\n');
end

% First, construct the transition matrix.
P = ConstructTransitionFromAdjacency(A);

% Parameter vectors, local gradient vectors, average gradient collections
G = zeros(d, n);
Z = zeros(d, n);
y = zeros(d,1);
fvals = zeros(T + 1, 1);
fits = zeros(n,1);
meanFits = zeros(T + 1,1);

for t = 1:T
  fvals(t) = ObjectiveFunction(y,optimal);
 %%
  for i=1:n
    fits(i) = ObjectiveFunction(X(:,i),optimal);
  end
  meanFits(t) = mean(fits);
 %%
 
  % Compute gradient for each objective function by taking inner product
  % with each column/row of B and X.  
 %%
 for i=1:n
     for j=1:d
        G(j,i) = GradientOfObjectiveFunction(X(:,i),optimal,j);
     end
 end
 %%%inprods = dot(B', X)';  
 %   G(:) = 0;
 %   G(:, (1 - inprods) > 0) = -B((1 - inprods) > 0, :)';
 %%
  
  Z = Z * P + G;
  % Project to 2-norm ball
  X = -stepsize * Z / sqrt(t);
  parameter_norms = sqrt(sum(X.^2));
  if (any(parameter_norms > R))
    large_indices = (parameter_norms > R);
    X(:, large_indices) = R * X(:, large_indices) ./ repmat(parameter_norms(large_indices), d, 1);
  end
  y = -stepsize * mean(Z, 2) / sqrt(t);
  if (norm(y) > R)
    y = R * y / norm(y);
  end
  if (mod(t, 100) == 0)
    fprintf(1, 'Distributed SVM Minimization Iteration %d, objval %2.3f\n', t, fvals(t));
  end
end

fvals(T + 1) = ObjectiveFunction(y,optimal);
%%
for i=1:n
    fits(i) = ObjectiveFunction(X(:,i),optimal);
end
meanFits(T+1) = mean(fits);
%%

