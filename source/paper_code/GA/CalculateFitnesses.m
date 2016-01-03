%Omid55
function [ fitnesses ] = CalculateFitnesses( population,o )

%% Fitness Caculation
%disp('CalculateFitnesses ... ');

fitnesses = zeros(size(population,1),1);
parfor i=1:size(population,1)
    fitnesses(i) = ObjectiveFunction(population(i,:),o);
end

end

