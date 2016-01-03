%Omid55
function [ population,fitnesses ] = Mutation4SVM( population,fitnesses,Pm,dataTrain,labelsTrain,dataTest,labelsTest )

disp('Mutation ... ');   
len = size(population,1);
t = 0;
mutatedPopulation = [];
for i=1:len
    
    t = t+1;
    index = randi(len,1);
    mutatedChromosome = population(index,:);
    
    r = rand(size(population,2),1);
    
    for j=1:size(r,1)
        if r(j) <= Pm
            mutatedChromosome(j) = mutatedChromosome(j) + 20 * rand(1,1) - 10;   %(rangeEnd - rangeBegin) * rand(1,1) - 1;
        end
    end
    mutatedPopulation = [mutatedPopulation; mutatedChromosome];

end

mutatedFitnesses = zeros(size(mutatedPopulation,1),1);
parfor i=1:size(mutatedPopulation,1)
    mutatedFitnesses(i) = CalculateFitness4SVM(mutatedPopulation,dataTrain,labelsTrain,dataTest,labelsTest);
end
population = [population; mutatedPopulation];
fitnesses = [fitnesses; mutatedFitnesses];


end

