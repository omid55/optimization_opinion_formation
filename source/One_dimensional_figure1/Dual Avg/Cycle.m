function A = Cycle(n, k)
% CYCLE Constructs adjacency matrix for n-cycle
%
% A = Cycle(n, k) sets A to be the symmetric adjacency matrix for the
%   k-connected n-cycle. Also adds self-loops. If k is not specified, sets
%   k = 1.

if (nargin < 2)
  k = 1;
end
if (k >= n)
  error('Must have n bigger than k');
end

A = eye(n);
for ii = 1:k
  A = A + circshift(eye(n), ii) + circshift(eye(n), -ii);
end
