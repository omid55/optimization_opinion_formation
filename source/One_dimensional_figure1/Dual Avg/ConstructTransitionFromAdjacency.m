function P = ConstructTransitionFromAdjacency(A)
% CONSTRUCTTRANSITIONFROMADJACENCY Constructs transition matrix from given
%   adjacency matrix.
%
% P = ConstructTransitionFromAdjacency(A) constructs the doubly stochastic
%   symmetric transition matrix from A by
%
%   P(i, j) = 1 / max(degree(i), degree(j))
%
%   and setting P(i, i) = 1 - sum_{i ~= j} P(i, j).

n = size(A, 1);
if (any(any(A - A' )))
  error('Adjacency matrix should be symmetric');
  return;
end
if (any(diag(A) == 0))
  for ii = 1:n
    if (A(ii, ii) == 0)
      A(ii, ii) = 1;
    end
  end
end

nneighbors = sum(A, 2);
Aleft = A ./ repmat(nneighbors, 1, n);
Aright = A ./ repmat(nneighbors', n, 1);
P = min(Aleft, Aright);
local_sums = sum(P) - diag(P)';
for ii = 1:n
  P(ii, ii) = 1 - local_sums(ii);
end
