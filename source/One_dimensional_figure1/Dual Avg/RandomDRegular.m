function A = RandomDRegular(d, n)
% RANDOMDREGULAR Constructs random d-regular graph
%
% A = RandomDRegular(d, n) constructs a random n-by-n adjacency matrix for a
%   graph that is d-regular. The current implementation is not completely
%   d-regular, but it is approximately d-regular.

numleft = d * ones(n, 1);
A = zeros(n);
for ii = 1:n
  possible_inds = find(numleft > 0);
  ind = find(possible_inds == ii);
  if (~isempty(ind))
    possible_inds = [possible_inds(1:ind-1) possible_inds(ind+1:end)];
  end
  selected = possible_inds(randperm(numel(possible_inds)));
  local_d = numleft(ii);
  selected = selected(1:min(local_d, numel(selected)));
  A(ii, selected) = 1;
  A(selected, ii) = 1;
  numleft(ii) = 0;
  numleft(selected) = numleft(selected) - 1;
end

A = A + eye(n);
