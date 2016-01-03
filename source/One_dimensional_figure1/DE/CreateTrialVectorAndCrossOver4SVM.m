%Omid55
function [ population ] = CreateTrialVectorAndCrossOver4SVM( population,fitnesses,Beta,Pr,Nv,dataTrain,labelsTrain,dataTest,labelsTest )

%disp('Create Trial Vector And Doing CrossOver ... ');
len = size(population,1);
nextGeneration = [];
for i=1:len
    
    Xi = population(i,:);
    fXi = fitnesses(i,:);
    
    diffs = zeros(1,size(population,2));
    for j=1:Nv
        list = 1:size(population,1);
        list(i) = [];
        xs = population(GetDistinctItems(list,3),:);
        diffs = diffs + xs(1,:) - xs(2,:);
    end
    Ui = xs(3,:) + Beta * diffs;
    
    % Binomial Crossover
    genesNum = size(population,2);
    jStar = randi(genesNum,1);
    J = jStar;
    for j=1:genesNum
        if rand() < Pr && j~=jStar
            J = [J j];
        end
    end
    
    child = Xi;
    child(J) = Ui(J);
    
    childFitness = CalculateFitness4SVM(child,dataTrain,labelsTrain,dataTest,labelsTest);
    if childFitness <= fXi
        nextGeneration = [nextGeneration; child];
    else
        nextGeneration = [nextGeneration; Xi];
    end

end

population = nextGeneration;


end
