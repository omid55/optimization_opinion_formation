%Omid55
%My Genetic Algorithm
%GA
function [  ] = GA( ranges,Use_Full_Crossover,Population_Model,id )

% init data
trainData=xlsread('\Data\Train.xlsx');  %train data => train inputs and outputs
trainData=trainData';
[R1,C1]=size(trainData);
trainInputs = trainData(1:R1-1,:);
trainTargets = trainData(R1,:);
testData=xlsread('\Data\Test.xlsx');  %test data => test inputs and outputs
testData=testData';
[R2,C2]=size(testData);
testInputs = testData(1:R2-1,:);
testTargets = testData(R2,:);

disp('--- GA Starts --- ');
%% Init Important Parameters
%Use_Full_Crossover = 0;           % 1 -> Full Crossover   and   0 -> Roullette Selection Crossover
%Population_Model = 0;               % 1 -> Steady State   and   0 -> Generational  
Classifier_Data_Limit = 500;
Number_Of_Population = 5;
Generation_Limit = 2;
Pc = 0.85;
Pm = 0.05;

% init variables
generation = 0;
generations = [];
s_generation = 0;
s_generations = [];
bests = [];
means = [];
s_bests = [];
s_means = [];
x = 0:0.01:1;


%% Body
population = InitPopulation(Number_Of_Population,ranges);
s_population = population;
fitnesses = CalculateFitnesses(population,trainInputs,trainTargets,testInputs,testTargets);
s_fitnesses = fitnesses;
learnerData = population;
learnerLabels = fitnesses;

% figure(1);
% populationV = GetPopulationValues(population);
% subplot(1,2,1);
% hist(populationV,x);

tic;
while TerminationSatisfy(s_fitnesses) ~= 1 && s_generation < Generation_Limit
    s_generation = s_generation + 1
    s_bests = [s_bests; max(s_fitnesses)];
    s_means = [s_means; mean(s_fitnesses)];
    s_generations = [s_generations; s_generation];

    [s_population,s_fitnesses,newData,newLabels] = CrossOver(0,s_population,s_fitnesses,Pc,Population_Model,trainInputs,trainTargets,testInputs,testTargets);

    [s_population,s_fitnesses,newData,newLabels] = Mutation(0,s_population,s_fitnesses,Pm,Population_Model,ranges,trainInputs,trainTargets,testInputs,testTargets);

    [s_population,s_fitnesses] = SelectionFunction(s_population,s_fitnesses,Number_Of_Population);
end
s_elapsedTime = toc;

tic;
while TerminationSatisfy(fitnesses) ~= 1 && generation < Generation_Limit
        generation = generation + 1
        bests = [bests; max(fitnesses)];
        means = [means; mean(fitnesses)];
        generations = [generations; generation];
    
        % Train The Learner
        disp('Training The Learner ...');
        if length(learnerLabels) > Classifier_Data_Limit
            learnerData = learnerData(1:Classifier_Data_Limit,:);
            learnerLabels = learnerLabels(1:Classifier_Data_Limit);
        end
        svmModel = svmtrain(learnerLabels, learnerData, '-s 3 -t 3 -r 1 -c 100 -b 1');
         
        if Use_Full_Crossover == 1
            % Full Crossover
            [population,fitnesses,newData,newLabels] = FullCrossOver(1,population,fitnesses,Pc,Population_Model,trainInputs,trainTargets,testInputs,testTargets,svmModel);
        else
            % Roulette Selection Crossover
            [s_population,s_fitnesses,newData,newLabels] = CrossOver(1,population,fitnesses,Pc,Population_Model,trainInputs,trainTargets,testInputs,testTargets,svmModel);
        end
        learnerData = [learnerData; newData];
        learnerLabels = [learnerLabels; newLabels];
        
        [population,fitnesses,newData,newLabels] = Mutation(1,population,fitnesses,Pm,Population_Model,ranges,trainInputs,trainTargets,testInputs,testTargets,svmModel);
        learnerData = [learnerData; newData];
        learnerLabels = [learnerLabels; newLabels];
        
        [population,fitnesses] = SelectionFunction(population,fitnesses,Number_Of_Population);
end
elapsedTime = toc;

% figure(1);
% populationV = GetPopulationValues(population);
% subplot(1,2,2);
% hist(populationV,x);


%% Ploting The Results
if Population_Model == 0
    popModel = 'Generational';
else
    popModel = 'Steady State';
end
if Use_Full_Crossover == 1
    crossoverType = 'Full';
else
    crossoverType = 'Roulette Selection';
end

figure;
leg = plot(generations,bests,'-v',generations,means,'-*',s_generations,s_bests,'-d',s_generations,s_means,'-p');
legend(leg,'Best of Our Algorithm Fitnesses','Mean of Our Algorithm Fitnesses','Best of Simple Fitnesses','Mean of Simple Fitnesses');
xlabel('Generations');
ylabel('Fitnesses');
title(sprintf('Simple Genetic in %d Iterations within %f Seconds & Our Genetic Algorithm in %d Iterations within %f Seconds with %s Population Model and %s Crossover by %d Chromosomes',length(s_generations),s_elapsedTime,length(generations),elapsedTime,popModel,crossoverType,Number_Of_Population));


save(sprintf('Data%d',id));


end

