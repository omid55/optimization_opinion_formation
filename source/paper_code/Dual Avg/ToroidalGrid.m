function A = ToroidalGrid(n, k, wrapover)
% TOROIDALGRID Constructs adjacency matrix of toroidal grid
%
% A = ToroidalGrid(n) returns the adjacency matrix of an n-by-n toroidal grid
%   in the matrix A. The matrix A has size n^2.
%
% A = ToroidalGrid(n, wrapover) is identical to the above, but if wrapover
%   is false, returns a non-toroidal grid.

if (nargin < 3)
  wrapover = true;
end

if(nargin < 2)
    k = 1;
end

if (wrapover)
  % Make adjacency matrix of k-connected cycle
  A = Cycle(n,k);
  % Take its graph cartesian product and add the diagonal. Note that cartesian
  % products of graphs are given by the kronecker sum of the adjacency
  % matrices.
  A = kron(A, eye(n)) + kron(eye(n), A) + eye(n^2);
else
  A = circshift(eye(n), 1) + circshift(eye(n), -1);
  A(1, n) = 0;
  A(n, 1) = 0;
  A = kron(A, eye(n)) + kron(eye(n), A) + eye(n^2);
end
