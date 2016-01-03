function [  ] = Run(  )
% TIMESCALINGEXPERIMENT Performs synthetic experiment to compare time
% scales of minimization
%
% TimeScalingExperiment compares the performance of the three distributed
% methods in SyntheticExperiments to one another in terms of time
% scaling. Fixes a graph of size 100 and compares cycle, grid, and expander
% in terms of time to get to a fixed level of accuracy.

close all;
clc;
%clear;

% Whether to do scaling by diameter or scaling by P's spectral gap
scaling_by_diameter = 0;

% Approxiate time to get to accuracy of 10^-2 for expander
T = 600;
% Number of nodes in graphs
n = 100;
% Number of tests per n
numtests = 3;
% Degree of random graph
degree = 5;
% Dimension of problems to solve
dim = 5;
% Radius of 2-norm ball
R = 5;
% Lipschitz constant
G = 1;
% Noise level to construct in synthetic data
noise_level = .1;

%% My Params
rangeBegin = -5;
rangeEnd = 5;
optimal = rangeBegin + (rangeEnd - rangeBegin) * rand(dim,1);
%%

A_grid = ToroidalGrid(round(sqrt(n)));
P_grid = ConstructTransitionFromAdjacency(A_grid);

A_cycle = Cycle(n);
P_cycle = ConstructTransitionFromAdjacency(A_cycle);

A_expander = RandomDRegular(degree, n);
P_expander = ConstructTransitionFromAdjacency(A_expander);

if (scaling_by_diameter)
  % The diameter of the random-d-regular graph seems to be about 5.
  T = 1000;
  T_grid = ceil(sqrt(n) / 5) * T;
  T_expander = T;
  T_cycle = ceil(n / 5) * T;
else
  % Scaling by spectral gap
  T_grid = ceil(T * (1 - norm(P_expander - ones(n) / n)) / (1 - norm(P_grid - ones(n) / n)));
  T_cycle = ceil(T * (1 - norm(P_expander - ones(n) / n)) / (1 - norm(P_cycle - ones(n) / n)) / 2);
  T_expander = T;
end

grid_opt_gaps = zeros(numtests, T_grid + 1);
cycle_opt_gaps = zeros(numtests, T_cycle + 1);
expander_opt_gaps = zeros(numtests, T_expander + 1);

eta = 4 * (R / G);

for test_iter = 1:numtests
  fprintf(1, 'Timing exp. with n = %d, iter %d of %d\n', n, test_iter, ...
          numtests);
  [opt_val B] = MakeSVMData(n, dim, noise_level, G, R);

  fprintf(1, '**** EXPANDER (%d iters) ****\n', T_expander);
  if (scaling_by_diameter)
    eta_expander = (1/5) * eta;
  else
    eta_expander = sqrt(1 - norm(P_expander - ones(n)/n)) * eta;
  end
  [obj_vals meanFitnesses] = DistributedSVM(B, A_expander, zeros(n, dim), T_expander, eta_expander, R, optimal);
  expander_opt_gaps(test_iter, :) = meanFitnesses;

%   fprintf(1, '**** CYCLE (%d iters) ****\n', T_cycle);
%   if (scaling_by_diameter)
%     eta_cycle = (1 / n) * eta;
%   else
%     eta_cycle = sqrt(1 - norm(P_cycle - ones(n)/n)) * eta;
%   end
%   [obj_vals meanFitnesses] = DistributedSVM(B, A_cycle, zeros(n,dim), T_cycle, eta_cycle, R, optimal);
%   cycle_opt_gaps(test_iter, :) = meanFitnesses;
end

figure;
plot(mean(expander_opt_gaps));
figure;
errorbar(mean(expander_opt_gaps),std(expander_opt_gaps));


end

