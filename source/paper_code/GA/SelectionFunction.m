%Omid55
function [ population,fitnesses ] = SelectionFunction( population,fitnesses,selectionLen )

sortedFitnesses = sort(fitnesses,'ascend');
threshold = sortedFitnesses(selectionLen);
population(fitnesses > threshold,:) = [];
fitnesses(fitnesses > threshold) = [];

indices = find(fitnesses == threshold);
removeList = GetDistinctItems(indices,size(population,1)-selectionLen);

removeList = sort(removeList,'descend');
population(removeList,:) = [];
fitnesses(removeList) = [];

end

