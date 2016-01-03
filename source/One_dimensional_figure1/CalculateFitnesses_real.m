%Omid55
%NewOrder-Based Crossovers for the Graph Coloring Problem
%Christine L. Mumford's Method

function [ fitnesses ] = CalculateFitnesses_real( population )

fitnesses = zeros(size(population,1),1);
for i=1:size(population,1)
    fitnesses(i) = CalculateFitness_real(population(i,:));
end

end

