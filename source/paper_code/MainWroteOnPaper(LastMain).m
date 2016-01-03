%Omid55
function [  ] = MainWroteOnPaper(  )

clc;
close all;

addpath('D:\Omid\matlab_bgl-4.0.1_2\matlab_bgl');
addpath('GA');
addpath('DE');
addpath('Dual Avg');

%% Parameters Initiallization
times = 2;
N = 100;
M = 5;
Mu = 0.9;
%Pa = 1;
MaxIteration = 2000;
%Eps = 10^-50;
rangeBegin = -5;
rangeEnd = 5;

ga_meanSum = zeros(times,MaxIteration);
de_meanSum = zeros(times,MaxIteration);
pso_meanSum = zeros(times,MaxIteration);
ba_meanSum = zeros(times,MaxIteration);
ws_meanSum = zeros(times,MaxIteration);
er_meanSum = zeros(times,MaxIteration);
dual_cycle_meanSum = zeros(times,MaxIteration);
dual_expander_meanSum = zeros(times,MaxIteration);

ga_bestSum = [];
de_bestSum = [];
pso_bestSum = [];
ba_bestSum = [];
ws_bestSum = [];
er_bestSum = [];
% ga_goodnessSum = [];
% de_goodnessSum = [];
% pso_goodnessSum = [];
% ba_goodnessSum = [];
% ws_goodnessSum = [];
% er_goodnessSum = [];

for time=1:times
    
    %tic;
    time

    optimal = rangeBegin + (rangeEnd - rangeBegin) * rand(1,M);
    result = repmat(optimal + zeros(1,M),N,1);        %  << CHECK HERE >>  FOR RESULT OF OBJECTIVE FUNCTION

    x = rangeBegin + (rangeEnd - rangeBegin) * rand(N,M);

%     P = randperm(size(x,2));
%     save('P_rand_Perm.mat','P');


   %% Dual Averaging Method  (DDA)
    % Whether to do scaling by diameter or scaling by P's spectral gap
    scaling_by_diameter = 0;
    % Approxiate time to get to accuracy of 10^-2 for expander
    T = 600;
    % Number of nodes in graphs
    %N = 100;
    % Degree of random graph
    degree = 5;
    % Dimension of problems to solve
    %M = 5;
    % Radius of 2-norm ball
    R = 5;
    % Lipschitz constant
    G = 1;
    % Noise level to construct in synthe%tic data
    noise_level = .1;

    % DDA - Cycle
    %A_cycle = Cycle(N);
    A_cycle = RandomDRegular(degree, N);
    A_expander = RandomDRegular(degree, N);
    P_cycle = ConstructTransitionFromAdjacency(A_cycle);
    P_expander = ConstructTransitionFromAdjacency(A_expander);
    
    if (scaling_by_diameter)
      % The diameter of the random-d-regular graph seems to be about 5.
      T = 1000;
      T_cycle = ceil(N / 5) * T;
    else
      % Scaling by spectral gap
      T_cycle = ceil(T * (1 - norm(P_expander - ones(N) / N)) / (1 - norm(P_cycle - ones(N) / N)) / 2);
    end

    cycle_opt_gaps = zeros(1, T_cycle + 1);

    eta = 4 * (R / G);

    [opt_val B] = MakeSVMData(N, M, noise_level, G, R);

    fprintf(1, '**** CYCLE (%d iters) ****\n', T_cycle);
    if (scaling_by_diameter)
      eta_cycle = (1 / N) * eta;
    else
      eta_cycle = sqrt(1 - norm(P_cycle - ones(N)/N)) * eta;
    end    
    %%%[obj_vals meanFits] = DistributedSVM(B, A_cycle, X, T_cycle, eta_cycle, R, optimal');
    [obj_vals meanFits] = DistributedSVM(B, A_cycle, x, MaxIteration-1, eta_cycle, R, optimal');
    
    meanFits = meanFits - 1;
    dual_results = obj_vals - 1;
    dual_cycle_meanSum(time,:) = meanFits;

    
   %% DDA Expander
    %A_expander = RandomDRegular(degree, N);
    A_expander = full(erdos_reyni(N,0.1));
    P_expander = ConstructTransitionFromAdjacency(A_expander);

    if (scaling_by_diameter)
      % The diameter of the random-d-regular graph seems to be about 5.
      T = 1000;
      T_expander = T;
    else
      % Scaling by spectral gap
      T_expander = T;
    end

    expander_opt_gaps = zeros(1, T_expander + 1);

    eta = 4 * (R / G);

    [opt_val B] = MakeSVMData(N, M, noise_level, G, R);

    fprintf(1, '**** EXPANDER (%d iters) ****\n', T_expander);
    if (scaling_by_diameter)
      eta_expander = (1/5) * eta;
    else
      eta_expander = sqrt(1 - norm(P_expander - ones(N)/N)) * eta;
    end
    %%%[obj_vals meanFits] = DistributedSVM(B, A_expander, X, T_expander, eta_expander, R, optimal');
    [obj_vals meanFits] = DistributedSVM(B, A_expander, x, MaxIteration-1, eta_expander, R, optimal');

    meanFits = meanFits - 1;
    dual_results = obj_vals - 1;
    dual_expander_meanSum(time,:) = meanFits;
    
    
   %% Simple GA
    disp('GA');
    Pc = 0.85;
    Pm = 0.05;
    ga_generation = 0;
    ga_population = x;
    ga_fitnesses = CalculateFitnesses(ga_population,optimal);
    ga_bests = [];
    ga_means = [];
    ga_generations = [];
    %ga_goodness = [];
    %ga_goodnessBest = [];
    %while TerminationSatisfy(ga_fitnesses,Eps) ~= 1 && ga_generation < MaxIteration
    while ga_generation < MaxIteration
        lastGa = ga_population;

        ga_generation = ga_generation + 1;
        ga_bests = [ga_bests; min(ga_fitnesses)];
        ga_means = [ga_means; mean(ga_fitnesses)];
        ga_generations = [ga_generations; ga_generation];
        %ga_goodness = [ga_goodness; norm(ga_population - result)];
        %bestIdx = find(ga_fitnesses == min(ga_fitnesses));
        %ga_goodnessBest = [ga_goodnessBest; norm(ga_population(bestIdx(1),:) - result(1,:))];

        [ga_population,ga_fitnesses] = CrossOver(ga_population,ga_fitnesses,Pc,optimal);

        [ga_population,ga_fitnesses] = Mutation(ga_population,ga_fitnesses,Pm,rangeBegin,rangeEnd,optimal);

        [ga_population,ga_fitnesses] = SelectionFunction(ga_population,ga_fitnesses,N);

%         if norm(ga_population-lastGa) < Eps
%             break;
%         end
    end
    
    ga_meanSum(time,:) = ga_means - 1;
    
    if size(ga_bestSum,1) == 0
        ga_bestSum = ga_bests;
    else
        ga_bestSum = ga_bestSum + ga_bests;
    end
%     if size(ga_goodnessSum,1) == 0
%         ga_goodnessSum = ga_goodness;
%     else
%         ga_goodnessSum = ga_goodnessSum + ga_goodness;
%     end


   %% PSO Inertia Weights
    disp('PSO');
    c1 = 0.1;
    c2 = 0.9;
    w = 0.001;
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
    %pso_goodness = [];
    %pso_goodnessBest = [];
    %while norm(pso_x-lastpso_x) > Eps && pso_it < MaxIteration
    while pso_it < MaxIteration

        lastpso_x = pso_x;
        pso_fitnesses = CalculateFitnesses(pso_x,optimal);
        pso_means = [pso_means; mean(pso_fitnesses)];
        pso_bests = [pso_bests; min(pso_fitnesses)];
        %pso_goodness = [pso_goodness; norm(pso_x - result)];
        %bestIdx = find(pso_fitnesses == min(pso_fitnesses));
        %pso_goodnessBest = [pso_goodnessBest; norm(pso_x(bestIdx(1),:) - result(1,:))];

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

    pso_meanSum(time,:) = pso_means - 1;
    
    if size(pso_bestSum,1) == 0
        pso_bestSum = pso_bests;
    else
        pso_bestSum = pso_bestSum + pso_bests;
    end
%     if size(pso_goodnessSum,1) == 0
%         pso_goodnessSum = pso_goodness;
%     else
%         pso_goodnessSum = pso_goodnessSum + pso_goodness;
%     end

    
   %% Differential Evolution (DE)
    disp('DE');
    Beta = 0.99;
    Pr = 0.99;
    Nv = 1;   % number of difference vectors
    % DE/rand/1/binomial
    de_population = x;
    fitnesses = CalculateFitnesses(de_population,optimal);
    de_it = 0;
    de_bests = [];
    de_means = [];
    %de_goodness = [];
    %de_goodnessBest = [];
    %while TerminationSatisfy(fitnesses,Eps) ~= 1 && de_it < MaxIteration
    while de_it < MaxIteration
        lastDe = de_population;

        de_it = de_it + 1;
        de_bests = [de_bests; min(fitnesses)];
        de_means = [de_means; mean(fitnesses)];
        %de_goodness = [de_goodness; norm(de_population-result)];
        %bestIdx = find(fitnesses == min(fitnesses));
        %de_goodnessBest = [de_goodnessBest; norm(de_population(bestIdx(1),:) - result(1,:))];

        [de_population] = CreateTrialVectorAndCrossOver(de_population,fitnesses,Beta,Pr,Nv,optimal);

        fitnesses = CalculateFitnesses(de_population,optimal);

%         if norm(de_population - lastDe) < Eps
%             break;
%         end
    end

    de_meanSum(time,:) = de_means - 1;
    
    if size(de_bestSum,1) == 0
        de_bestSum = de_bests;
    else
        de_bestSum = de_bestSum + de_bests;
    end
%     if size(de_goodnessSum,1) == 0
%         de_goodnessSum = de_goodness;
%     else
%         de_goodnessSum = de_goodnessSum + de_goodness;
%     end    


   %% My Algorithm's Simulation Body (Consensus Based Optimization)(CBO) 
    % CBO BA 
    disp('CBO-BA');

   %% Network creation
    %net = BarabasiGraphCreator(N,4,4);
    
    averageDegree = 12;
    net = BarabasiGraphCreator(N,averageDegree/2,averageDegree/2);
   %%
   
    originalX = x;
    
    %ba_diffs = [];
    ba_bests = [];
    ba_means = [];
    %ba_myGoodness = [];
    %ba_myGoodnessBest = [];
    fits = zeros(size(x,1),1);
    ba_xs = zeros(MaxIteration*N,M);
    for l=1:MaxIteration
        % just for ploting the graph
        for i=1:size(x,1)
            fits(i) = ObjectiveFunction(x(i,:),optimal);
        end
        ba_means = [ba_means; mean(fits)];
        ba_bests = [ba_bests; min(fits)];
%        bestIdx = find(fits == min(fits));
%         ba_myGoodnessBest = [ba_myGoodnessBest; norm(x(bestIdx(1),:) - result(1,:))];
%         ba_myGoodness = [ba_myGoodness; norm(x - result)];
        % just for plotting the graph

        lastX = x;

        for i=1:N  % for each agent
            
            adj = net(num2str(i));
            if length(adj) == 0
                break;
            end

            dummy = 0.2;      % << CHECK HERE >>   FOR VALUE OF DUMMY VARIABLE
            for feat=1:M
                fitns = zeros(length(adj),1);
                for j=1:length(adj)
                    mask = dummy * ones(1,M);   % dummy values
                    mask(feat) = x(adj(j),feat);
                    fitns(j) = ObjectiveFunction(mask,optimal);
                end
                bestIdx = find(fitns == min(fitns));
                bestAdj = adj(bestIdx(1));
                selfMask = dummy * ones(1,M);
                selfMask(feat) = x(i,feat);
                selfMaskedFit = ObjectiveFunction(selfMask,optimal);
                %w = 1 / min(fits);
% % %                 w = selfMaskedFit / min(fits);
% % %                 if selfMaskedFit > min(fits) || Pa < rand(1)    % Pa==1 -> if exists   and   Pa==0 -> if doesn't exist
                w = min(fitns) / selfMaskedFit;
                
                if w<1
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

%         diff = norm(x - lastX);
%         ba_diffs = [ba_diffs; diff];
%         if diff < Eps
%             break;
%         end

        ba_xs((l-1)*N+1:l*N,:) = x;
    end
    % End Of CBO-BA
    
    ba_meanSum(time,:) = ba_means - 1;
    
    if size(ba_bestSum,1) == 0
        ba_bestSum = ba_bests;
    else
        ba_bestSum = ba_bestSum + ba_bests;
    end
%    if size(ba_goodnessSum,1) == 0
%         ba_goodnessSum = ba_myGoodness;
%     else
%         ba_goodnessSum = ba_goodnessSum + ba_myGoodness;
%     end

    
   %% CBO WS
    disp('CBO-WS');

   %% Network creation
    %sp = WattsStrogatzCreator(N,4,0.1);
    
    sp = WattsStrogatzCreator(N,averageDegree,0.1);
    net = CreateMap(sp);
   %%
    
    x = originalX;
    
    %ws_diffs = [];
    ws_bests = [];
    ws_means = [];
    %ws_myGoodness = [];
    %ws_myGoodnessBest = [];
    ws_xs = zeros(MaxIteration*N,M);
    fits = zeros(size(x,1),1);
    for l=1:MaxIteration

        % just for ploting the graph
        for i=1:size(x,1)
            fits(i) = ObjectiveFunction(x(i,:),optimal);
        end
        ws_means = [ws_means; mean(fits)];
        ws_bests = [ws_bests; min(fits)];
%        bestIdx = find(fits == min(fits));
%         ws_myGoodnessBest = [ws_myGoodnessBest; norm(x(bestIdx(1),:) - result(1,:))];
%         ws_myGoodness = [ws_myGoodness; norm(x - result)];
        % just for plotting the graph

        lastX = x;

        for i=1:N  % for each agent
            
            adj = net(num2str(i));
            if length(adj) == 0
                break;
            end

            dummy = 0.2;      % << CHECK HERE >>   FOR VALUE OF DUMMY VARIABLE
            for feat=1:M
                fitns = zeros(length(adj),1);
                for j=1:length(adj)
                    mask = dummy * ones(1,M);   % dummy values
                    mask(feat) = x(adj(j),feat);
                    fitns(j) = ObjectiveFunction(mask,optimal);
                end
                bestIdx = find(fitns == min(fitns));
                bestAdj = adj(bestIdx(1));
                selfMask = dummy * ones(1,M);
                selfMask(feat) = x(i,feat);
                selfMaskedFit = ObjectiveFunction(selfMask,optimal);
                %w = 1 / min(fits);
% % %                 w = selfMaskedFit / min(fits);
% % %                 if selfMaskedFit > min(fits) || Pa < rand(1)    % Pa==1 -> if exists   and   Pa==0 -> if doesn't exist
                w = min(fitns) / selfMaskedFit;
                
                if w<1
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

%         diff = norm(x - lastX);
%         ws_diffs = [ws_diffs; diff];
%         if diff < Eps
%             break;
%         end

        ws_xs((l-1)*N+1:l*N,:) = x;
    end
    % End Of CBO-WS
    
    ws_meanSum(time,:) = ws_means - 1;
    
    if size(ws_bestSum,1) == 0
        ws_bestSum = ws_bests;
    else
        ws_bestSum = ws_bestSum + ws_bests;
    end
    

%     if size(ws_goodnessSum,1) == 0
%         ws_goodnessSum = ws_myGoodness;
%     else
%         ws_goodnessSum = ws_goodnessSum + ws_myGoodness;
%     end
    
    
   %% CBO ER
    disp('CBO-ER');

   %% Network creation
    %erdosProbability = 0.1;
   
    erdosProbability = (50+averageDegree)/N;   % more average degree for ER makes better solutions
    sp = erdos_reyni(N,erdosProbability);
    net = CreateMap(sp);
   %%
    
    x = originalX;
    
    %er_diffs = [];
    er_bests = [];
    er_means = [];
    %er_myGoodness = [];
    %er_myGoodnessBest = [];
    er_xs = zeros(MaxIteration*N,M);
    fits = zeros(size(x,1),1);
    for l=1:MaxIteration

        % just for ploting the graph
        for i=1:size(x,1)
            fits(i) = ObjectiveFunction(x(i,:),optimal);
        end
        er_means = [er_means; mean(fits)];
        er_bests = [er_bests; min(fits)];
%        bestIdx = find(fits == min(fits));
%         er_myGoodnessBest = [er_myGoodnessBest; norm(x(bestIdx(1),:) - result(1,:))];
%         er_myGoodness = [er_myGoodness; norm(x - result)];
        % just for plotting the graph

        lastX = x;

        for i=1:N  % for each agent
            
            adj = net(num2str(i));
            if length(adj) == 0
                break;
            end

            dummy = 0.2;      % << CHECK HERE >>   FOR VALUE OF DUMMY VARIABLE
            for feat=1:M
                fitns = zeros(length(adj),1);
                for j=1:length(adj)
                    mask = dummy * ones(1,M);   % dummy values
                    mask(feat) = x(adj(j),feat);
                    fitns(j) = ObjectiveFunction(mask,optimal);
                end
                bestIdx = find(fitns == min(fitns));
                bestAdj = adj(bestIdx(1));
                selfMask = dummy * ones(1,M);
                selfMask(feat) = x(i,feat);
                selfMaskedFit = ObjectiveFunction(selfMask,optimal);
                %w = 1 / min(fits);
% % %                 w = selfMaskedFit / min(fits);
% % %                 if selfMaskedFit > min(fits) || Pa < rand(1)    % Pa==1 -> if exists   and   Pa==0 -> if doesn't exist
                w = min(fitns) / selfMaskedFit;
                
                if w<1
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

%         diff = norm(x - lastX);
%         er_diffs = [er_diffs; diff];
%         if diff < Eps
%             break;
%         end

        er_xs((l-1)*N+1:l*N,:) = x;
    end
    % End Of CBO-ER
    
    er_meanSum(time,:) = er_means - 1;
    
    if size(er_bestSum,1) == 0
        er_bestSum = er_bests;
    else
        er_bestSum = er_bestSum + er_bests;
    end
%    if size(er_goodnessSum,1) == 0
%         er_goodnessSum = er_myGoodness;
%     else
%         er_goodnessSum = er_goodnessSum + er_myGoodness;
%     end

    %toc
    ooo = 0;
    
end

% ga_bestSum = ga_bestSum / times;
% de_bestSum = de_bestSum / times;
% pso_bestSum = pso_bestSum / times;
% ba_bestSum = ba_bestSum / times;
% ws_bestSum = ws_bestSum / times;
% er_bestSum = er_bestSum / times;
% ga_goodnessSum = ga_goodnessSum / times;
% de_goodnessSum = de_goodnessSum / times;
% pso_goodnessSum = pso_goodnessSum / times;
% ba_goodnessSum = ba_goodnessSum / times;
% ws_goodnessSum = ws_goodnessSum / times;
% er_goodnessSum = er_goodnessSum / times;

FName = 'F2';
fig = figure;
hold on;
set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
errorbar(1:MaxIteration,mean(dual_cycle_meanSum),std(dual_cycle_meanSum),'c');
errorbar(1:MaxIteration,mean(dual_expander_meanSum),std(dual_expander_meanSum),'c');
errorbar(1:MaxIteration,mean(ga_meanSum),std(ga_meanSum),'r');
errorbar(1:MaxIteration,mean(de_meanSum),std(de_meanSum),'b');
errorbar(1:MaxIteration,mean(pso_meanSum),std(pso_meanSum),'g');
errorbar(1:MaxIteration,mean(ba_meanSum),std(ba_meanSum),'m');
errorbar(1:MaxIteration,mean(ws_meanSum),std(ws_meanSum),'y');
errorbar(1:MaxIteration,mean(er_meanSum),std(er_meanSum),'k');
leg = legend('DDA-C','DDA-E','GA','DE','PSO','CBO-BA','CBO-WS','CBO-ER','Location','Northeast');
set(leg,'FontSize',7);
xlabel('\bfIteration #');
ylabel(['\bf',FName]);
set(gca,'XScale','log');
%set(gca,'YScale','log');
saveas(fig,'MeansErrorBar.fig');
print -loos -dtiff MeansErrorBar.tiff;

%FName = 'F2';
fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
plot(1:MaxIteration,mean(dual_cycle_meanSum),'-o',1:MaxIteration,mean(dual_expander_meanSum),'-^',1:MaxIteration,mean(ga_meanSum),'-s',1:MaxIteration,mean(de_meanSum),'-d',1:MaxIteration,mean(pso_meanSum),'-x',1:MaxIteration,mean(ba_meanSum),'-p',1:MaxIteration,mean(ws_meanSum),'-v',1:MaxIteration,mean(er_meanSum),'-*');
leg = legend('DDA-C','DDA-E','GA','DE','PSO','CBO-BA','CBO-WS','CBO-ER','Location','Northeast');
set(leg,'FontSize',7);
xlabel('\bfIteration #');
ylabel(['\bf',FName]);
set(gca,'XScale','log');
%set(gca,'YScale','log');
saveas(fig,'Means.fig');
print -loos -dtiff Means.tiff;


save('F2DataFinal.mat');



end

