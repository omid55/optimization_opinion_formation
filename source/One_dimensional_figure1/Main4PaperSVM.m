%Omid55
%Optimization With Opinion Formation Models
function [  ] = Main4PaperSVM(  )


%% Clear All & Close All
clc;
close all;


%% Initiallization
N = 5;
M = 4;
Mu = 0.9;
Pa = 1;
MaxIteration = 3;
Eps = 10^-50;
rangeBegin = 0;
rangeEnd = 100;

% % for nn training
% rangeBegin = 0.1;
% rangeEnd = 1;


%% Creating Network
isClique = 0;

%NETS

netModel = 'BA';
net = BarabasiGraphCreator(N,5,5);

% netModel = 'ER';
% averageDegree = 30;
% erdosProbability = averageDegree/N;
% sp = erdos_reyni(N,erdosProbability);
% mat = CreateMap(sp);

% netModel = 'WS';
% sp = WMStrogatzCreator(N,4,0.1);
% mat = CreateMap(sp);

% netModel = 'Regu;ar';
% sp = CreateRegularLattice(N,4);
% mat = CreateMap(sp);
 
% netModel = 'Clique';
% isClique = 1;


%% Calculating Values
x = rangeBegin + (rangeEnd - rangeBegin) * rand(N,M);


%% Data Loading
addpath('libsvm');
load iris.dat;
data = iris(randperm(size(iris)),:);
trainLen = ceil(0.7 * size(data,1));
dataTrain = data(1:trainLen,1:end-1);
labelsTrain = data(1:trainLen,end);
dataTest = data(trainLen+1:end,1:end-1);
labelsTest = data(trainLen+1:end,end);
save('MyData.mat','dataTrain','labelsTrain','dataTest','labelsTest');


%% Simple GA
addpath('GA');
tic;
Pc = 0.85;
Pm = 0.05;
ga_generation = 0;
ga_population = x;
ga_fitnesses = zeros(size(ga_population,1),1);
parfor i=1:size(ga_population,1)
    ga_fitnesses(i) = ObjectiveFunction4SVM(ga_population(i,:),dataTrain,labelsTrain,dataTest,labelsTest);
end
ga_bests = [];
ga_means = [];
ga_generations = [];
ga_goodness = [];
ga_goodnessBest = [];
%while TerminationSatisfy(ga_fitnesses,Eps) ~= 1 && ga_generation < MaxIteration
while ga_generation < MaxIteration
    lastGa = ga_population;
    
    ga_generation = ga_generation + 1
    ga_bests = [ga_bests; max(ga_fitnesses)];
    ga_means = [ga_means; mean(ga_fitnesses)];
    ga_generations = [ga_generations; ga_generation];
    bestIdx = find(ga_fitnesses == min(ga_fitnesses));

    [ga_population,ga_fitnesses] = CrossOver4SVM(ga_population,ga_fitnesses,Pc,dataTrain,labelsTrain,dataTest,labelsTest);

    [ga_population,ga_fitnesses] = Mutation4SVM(ga_population,ga_fitnesses,Pm,dataTrain,labelsTrain,dataTest,labelsTest);

    [ga_population,ga_fitnesses] = SelectionFunction(ga_population,ga_fitnesses,N);
    
    if norm(ga_population-lastGa) < Eps
        break;
    end
end
ga_accuracy = 100 - 100 * ga_fitnesses(bestIdx(1),:);
ga_elapsedTime = toc;


%% PSO Inertia Weights
tic;
c1 = 0.9;
c2 = 0.1;
w = 0.1;
%w=1;                    << Original PSO >>

pso_x = x;
y = zeros(N,M);
yCosts = 99999 * ones(N,1);
yHat = zeros(1,M);
yHatCost = 99999;
v = zeros(N,M);
lastpso_x = zeros(N,M);

pso_it = 0;
pso_means = [];
pso_bests = [];
pso_goodnessBest = [];
while norm(pso_x-lastpso_x) > Eps && pso_it < MaxIteration
    if mod(pso_it,100)  == 0
        pso_it
    end

    lastpso_x = pso_x;
    pso_fitnesses = zeros(size(pso_x,1),1);
    parfor i=1:size(pso_x,1)
        pso_fitnesses(i) = ObjectiveFunction4SVM(pso_x(i,:),dataTrain,labelsTrain,dataTest,labelsTest);
    end
    
    pso_means = [pso_means; mean(pso_fitnesses)];
    pso_bests = [pso_bests; min(pso_fitnesses)];
    bestIdx = find(pso_fitnesses == min(pso_fitnesses));
    
    % Yi and YHat Calculation
    for i=1:N
        pso_xCost = pso_fitnesses(i);
        if pso_xCost < yCosts(i)
            y(i,:) = pso_x(i,:);
            yCosts(i) = pso_xCost;
        end
        if yCosts(i) < yHatCost
            yHat = y(i,:);
            yHatCost = yCosts(i);
        end
    end

    % updating Vi and pso_xi
    for i=1:N
        v(i,:) = w * v(i,:) + c1 * rand(1,M) .* (y(i,:) - pso_x(i,:)) + c2 * rand(1,M) .* (yHat - pso_x(i,:));
    end
    pso_x = pso_x + v;
    
    pso_it = pso_it + 1;
    
end
pso_accuracy = 100 - 100 * pso_fitnesses(bestIdx(1),:);
pso_elapsedTime = toc;


%% Differential Evolution (DE)
disp('DE');
addpath('DE');
tic;
Beta = 0.5;
Pr = 0.3;
Nv = 1;   % number of difference vectors
% DE/rand/1/binomial
de_population = x;
de_fitnesses = zeros(size(de_population,1),1);
parfor i=1:size(de_population,1)
    de_fitnesses(i) = ObjectiveFunction4SVM(de_population(i,:),dataTrain,labelsTrain,dataTest,labelsTest);
end
de_it = 0;
de_bests = [];
de_means = [];
de_goodnessBest = [];
%while TerminationSatisfy(fitnesses,Eps) ~= 1 && de_it < MaxIteration
while de_it < MaxIteration
    lastDe = de_population;
    
    de_it = de_it + 1;
    de_bests = [de_bests; min(de_fitnesses)];
    de_means = [de_means; mean(de_fitnesses)];
    bestIdx = find(de_fitnesses == min(de_fitnesses));
    
    [de_population] = CreateTrialVectorAndCrossOver4SVM(de_population,de_fitnesses,Beta,Pr,Nv,dataTrain,labelsTrain,dataTest,labelsTest);
 
    parfor i=1:size(de_population,1)
        de_fitnesses(i) = ObjectiveFunction4SVM(de_population(i,:),dataTrain,labelsTrain,dataTest,labelsTest);
    end
    
    if norm(de_population - lastDe) < Eps
        break;
    end
end
de_accuracy = 100 - 100 * de_fitnesses(bestIdx(1),:);
de_elapsedTime = toc;


%% My Algorithm's Simulation Body (Consensus Based Optimization)(CBO)
disp('CBO');

tic;
diffs = [];
bests = [];
means = [];
xs = [];
for l=1:MaxIteration
    if mod(l,10)==0
        l
    end
    
    % just for ploting the graph
    fits = zeros(size(x,1),1);
    for i=1:size(x,1)
        fits(i) = ObjectiveFunction4SVM(x(i,:),dataTrain,labelsTrain,dataTest,labelsTest);
    end
    means = [means; mean(fits)];
    bests = [bests; min(fits)];
    bestIdx = find(fits == min(fits));
    % just for plotting the graph
    
    lastX = x;
    
    for i=1:N  % for each agent
        adj = net(num2str(i));
        if length(adj) == 0
            break;
        end
                
        dummy = 0.2;      % << CHECK HERE >>   FOR VALUE OF DUMMY VARIABLE
        for feat=1:M
            fits = zeros(length(adj),1);
            for j=1:length(adj)
                mask = dummy * ones(1,M);   % dummy values
                mask(feat) = x(adj(j),feat);
                fits(j) = ObjectiveFunction4SVM(mask,dataTrain,labelsTrain,dataTest,labelsTest);
            end
            bestIdx = find(fits == min(fits));
            bestAdj = adj(bestIdx(1));
            selfMask = dummy * ones(1,M);
            selfMask(feat) = x(i,feat);
            selfMaskedFit = ObjectiveFunction4SVM(selfMask,dataTrain,labelsTrain,dataTest,labelsTest);
            %w = 1 / min(fits);
            w = selfMaskedFit / min(fits);
            if selfMaskedFit > min(fits) || Pa < rand(1)    % Pa==1 -> if exists   and   Pa==0 -> if doesn't exist
                x(i,feat) = x(i,feat) + Mu * w * (x(bestAdj,feat) - x(i,feat));
                if x(i,feat) > rangeEnd
                    x(i,feat) = rangeEnd;
                else 
                    if x(i,feat) < rangeBegin
                        x(i,feat) = rangeBegin;
                    end
                end
            end
        end
        
    end
    
    diff = norm(x - lastX);
    diffs = [diffs; diff];
    if diff < Eps
        break;
    end
    
    xs = [xs; x];
end

LastIteration = l
cbo_accuracy =  100 - 100 * fits(bestIdx(1),:);


% ------------ === FIGURES === ------------
%------------------------------------------------------------------------------------------------------------
f = figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
for i=1:size(xs,2)
    subplot(1,2,1),plot(xs(:,i));
    hold on;
end
%title('\bfOpinion Changing');
title('\bf(a)');
xlabel('\bfIteration #');
ylabel('\bfOpinion');

subplot(1,2,2),plot(diffs);
%title('\bfGrow of Opinions');
title('\bf(b)');
xlabel('\bfIteration #');
ylabel('\bfDifference');
saveas(f,'OpinionsChanges.fig');
print -loos -dtiff OpinionsChanges.tiff;
%------------------------------------------------------------------------------------------------------------

%=============================================================================
figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
plot(1:LastIteration,means,'-*',1:LastIteration,bests,'-x');
legend('mean','best');
xlabel('\bfIteration #');
ylabel('\bfValue');
saveas(f,'meanbest.fig');
print -loos -dtiff meanbest.tiff;
%=============================================================================


%------------------------------------------------------------------------------------------------------------
fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
plot(ga_generations,ga_means,'-*',1:de_it,de_means,'-^',1:pso_it,pso_means,'-x',1:LastIteration,means,'-p');
legend('GA','DE','PSO','CBO','Location','best');
xlabel('\bfIteration #');
ylabel('\bfValue');
title('\bfMean of Objective Function');

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
plot(ga_generations,ga_bests,'-*',1:de_it,de_bests,'-^',1:pso_it,pso_bests,'-x',1:LastIteration,bests,'-p');
legend('GA','DE','PSO','CBO','Location','best');
xlabel('\bfIteration #');
ylabel('\bfValue');
title('\bfBest of Objective Function');
%------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------------------
fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
plot(ga_generations,log(ga_means),'-*',1:de_it,log(de_means),'-^',1:pso_it,log(pso_means),'-x',1:LastIteration,log(means),'-p');
legend('GA','DE','PSO','CBO','Location','best');
xlabel('\bfIteration #');
ylabel('\bfLog Value');
title('\bfLog Mean of Objective Function');

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
plot(ga_generations,log(ga_bests),'-*',1:de_it,log(de_bests),'-^',1:pso_it,log(pso_bests),'-x',1:LastIteration,log(bests),'-p');
legend('GA','DE','PSO','CBO','Location','best');
xlabel('\bfIteration #');
ylabel('\bfValue');
title('\bfLog Best of Objective Function');
%------------------------------------------------------------------------------------------------------------

ga_accuracy
de_accuracy
pso_accuracy
cbo_accuracy


save('TheWholeData.mat');
% pause(30);
% exit;


end

