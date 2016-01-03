%Omid55
function [  ] = Main(  )

clc;
close all;

addpath('GA');
addpath('DE');

times = 50;
N = 100;
M = 1;
Mu = 0.9;
MaxIteration = 50;
Eps = 10^-50;
rangeBegin = -1;
rangeEnd = 1;

ga_meanSum = zeros(times,MaxIteration);
de_meanSum = zeros(times,MaxIteration);
pso_meanSum = zeros(times,MaxIteration);
ba_meanSum = zeros(times,MaxIteration);
ws_meanSum = zeros(times,MaxIteration);
er_meanSum = zeros(times,MaxIteration);
dual_cycle_meanSum = zeros(times,MaxIteration);
dual_expander_meanSum = zeros(times,MaxIteration);

for time=1:times
    
    time
    
    % % REAL NETWORK
    % network_file = 'USairport500-net.txt';
    % a = importdata(network_file);
    % MM = max(max(a(:,1)),max(a(:,2)));
    % sp = sparse(a(:,1),a(:,2),1,MM,MM);
    % net = CreateMap(sp);
    % N = size(sp,1);
    % % REAL NETWORK

    x = rangeBegin + (rangeEnd - rangeBegin) * rand(N,1);
    x = InitPopulation_real(x);
    M = size(x,2);

    % P = randperm(size(x,2));
    % save('P_rand_Perm.mat','P');
    

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
    A_cycle = Cycle(N);
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
    %%%[obj_vals meanFits] = DistributedSVM(B, A_cycle, X, T_cycle, eta_cycle, R);
    [obj_vals meanFits] = DistributedSVM(B, A_cycle, x, MaxIteration-1, eta_cycle, R);
    
    dual_cycle_meanSum(time,:) = meanFits;

    
   %% DDA Expander
    A_expander = RandomDRegular(degree, N);
    %A_expander = full(erdos_reyni(N,0.1));
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
    %%%[obj_vals meanFits] = DistributedSVM(B, A_expander, X, T_expander, eta_expander, R);
    [obj_vals meanFits] = DistributedSVM(B, A_expander, x, MaxIteration-1, eta_expander, R);

    dual_expander_meanSum(time,:) = meanFits;
    

   %% Simple GA
    disp('GA');
    tic;
    Pc = 0.85;
    Pm = 0.05;
    ga_generation = 0;
    ga_population = x;
    ga_fitnesses = CalculateFitnesses_real(ga_population);
    ga_bests = [];
    ga_means = [];
    ga_generations = [];
    %while TerminationSatisfy(ga_fitnesses,Eps) ~= 1 && ga_generation < MaxIteration
    while ga_generation < MaxIteration
        lastGa = ga_population;

        ga_generation = ga_generation + 1;
        ga_bests = [ga_bests; max(ga_fitnesses)];
        ga_means = [ga_means; mean(ga_fitnesses)];
        
        ga_generations = [ga_generations; ga_generation];
        bestIdx = find(ga_fitnesses == min(ga_fitnesses));

        [ga_population,ga_fitnesses] = CrossOver(ga_population,ga_fitnesses,Pc);

        [ga_population,ga_fitnesses] = Mutation(ga_population,ga_fitnesses,Pm,rangeBegin,rangeEnd);

        [ga_population,ga_fitnesses] = SelectionFunction(ga_population,ga_fitnesses,N);

%         if norm(ga_population-lastGa) < Eps
%             break;
%         end
    end
    ga_elapsedTime = toc;
    
    ga_meanSum(time,:) = ga_means - 0.5;
    

   %% PSO Inertia Weights
    disp('PSO');
    tic;
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
    while pso_it < MaxIteration
    %while norm(pso_x-lastpso_x) > Eps && pso_it < MaxIteration

        lastpso_x = pso_x;
        pso_fitnesses = CalculateFitnesses_real(pso_x);
        pso_means = [pso_means; mean(pso_fitnesses)];
        pso_bests = [pso_bests; min(pso_fitnesses)];

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
    pso_elapsedTime = toc;
    
    pso_meanSum(time,:) = pso_means - 0.5;
    

   %% Differential Evolution (DE)
    disp('DE');
    tic;
    Beta = 0.99;
    Pr = 0.99;
    Nv = 1;   % number of difference vectors
    % DE/rand/1/binomial
    de_population = x;
    fitnesses = CalculateFitnesses_real(de_population);
    de_it = 0;
    de_bests = [];
    de_means = [];
    %while TerminationSatisfy(fitnesses,Eps) ~= 1 && de_it < MaxIteration
    while de_it < MaxIteration
        lastDe = de_population;

        de_it = de_it + 1;
        de_bests = [de_bests; min(fitnesses)];
        de_means = [de_means; mean(fitnesses)];
        [de_population] = CreateTrialVectorAndCrossOver(de_population,fitnesses,Beta,Pr,Nv);

        fitnesses = CalculateFitnesses_real(de_population);

%         if norm(de_population - lastDe) < Eps
%             break;
%         end
    end
    de_elapsedTime = toc;

    de_meanSum(time,:) = de_means - 0.5;


  %% My Algorithm's Simulation Body (Consensus Based Optimization)(CBO)
   originalX = x;
  
   % CBO BA 
   disp('CBO-BA');
   
   %% Network creation
    %net = BarabasiGraphCreator(N,4,4);
    
    averageDegree = 12;
    net = BarabasiGraphCreator(N,averageDegree/2,averageDegree/2);
   %%
   
    %ba_bests = [];
    ba_means = [];
    %ba_xs = [];
    for l=1:MaxIteration

        % just for ploting the graph
        fits = zeros(size(x,1),1);
        for i=1:size(x,1)
            fits(i) = CalculateFitnesses_real(x(i,:));
        end
        ba_means = [ba_means; mean(fits)];
        %ba_bests = [ba_bests; min(fits)];
        % just for plotting the graph

        lastX = x;

        for i=1:N  % for each agent
            
            adj = net(num2str(i));
            if length(adj) == 0
                break;
            end

            adjFits = zeros(length(adj),1);
            for j=1:length(adj)
                adjFits(j) = CalculateFitnesses_real(x(adj(j),:));
            end
            idx = find(adjFits == min(adjFits));
            %selfF = CalculateFitnesses_real(x(i,:));
            
            w = 0.5/min(adjFits);
            %w = selfF / min(adjFits);
            bestPal = adj(idx);
            x(i,:) = x(i,:) + Mu * w * (x(bestPal(1),:) - x(i,:));
            x(i,x(i,:) > 1) = 1;
            x(i,x(i,:) < -1) = -1;               
            
        end

%         diff = norm(x - lastX);
%         ba_diffs = [ba_diffs; diff];
%         if diff < Eps
%             break;
%         end

        %ba_xs = [ba_xs; x];
    end
    % End Of CBO-BA

    ba_meanSum(time,:) = ba_means - 0.5;
    
    
   %% CBO WS
    disp('CBO-WS');

   %% Network creation
    %sp = WattsStrogatzCreator(N,4,0.1);
    
    sp = WattsStrogatzCreator(N,averageDegree,0.1);
    net = CreateMap(sp);
   %%    
    
    x = originalX;
    
    ws_diffs = [];
    ws_bests = [];
    ws_means = [];
    ws_xs = [];
    for l=1:MaxIteration

        % just for ploting the graph
        fits = zeros(size(x,1),1);
        for i=1:size(x,1)
            fits(i) = CalculateFitnesses_real(x(i,:));
        end
        ws_means = [ws_means; mean(fits)];
        ws_bests = [ws_bests; min(fits)];
        bestIdx = find(fits == min(fits));
        % just for plotting the graph

        lastX = x;

        for i=1:N  % for each agent
            
            adj = net(num2str(i));
            if length(adj) == 0
                break;
            end

            adjFits = zeros(length(adj),1);
            for j=1:length(adj)
                adjFits(j) = CalculateFitnesses_real(x(adj(j),:));
            end
            idx = find(adjFits == min(adjFits));
            %selfF = CalculateFitnesses_real(x(i,:));
            
            w = 0.5/min(adjFits);
           %w = selfF / min(adjFits);
            bestPal = adj(idx);
            x(i,:) = x(i,:) + Mu * w * (x(bestPal(1),:) - x(i,:));
            x(i,x(i,:) > 1) = 1;
            x(i,x(i,:) < -1) = -1;
            
        end

        diff = norm(x - lastX);
        ws_diffs = [ws_diffs; diff];
%         if diff < Eps
%             break;
%         end

        ws_xs = [ws_xs; x];
    end
    % End Of CBO-WS

    ws_meanSum(time,:) = ws_means - 0.5;
    
    
   %% CBO ER
    disp('CBO-ER');

   %% Network creation
    %erdosProbability = 0.1;
   
    erdosProbability = (30+averageDegree)/N;   % more average degree for ER makes better solutions
    sp = erdos_reyni(N,erdosProbability);
    net = CreateMap(sp);
   %%
    
    x = originalX;
    
    er_diffs = [];
    er_bests = [];
    er_means = [];
    er_xs = [];
    for l=1:MaxIteration

        % just for ploting the graph
        fits = zeros(size(x,1),1);
        for i=1:size(x,1)
            fits(i) = CalculateFitnesses_real(x(i,:));
        end
        er_means = [er_means; mean(fits)];
        er_bests = [er_bests; min(fits)];
        bestIdx = find(fits == min(fits));
        % just for plotting the graph

        lastX = x;

        for i=1:N  % for each agent
            
            adj = net(num2str(i));
            if length(adj) == 0
                break;
            end

            adjFits = zeros(length(adj),1);
            for j=1:length(adj)
                adjFits(j) = CalculateFitnesses_real(x(adj(j),:));
            end
            idx = find(adjFits == min(adjFits));
            %selfF = CalculateFitnesses_real(x(i,:));
            
            w = 0.5/min(adjFits);
           %w = selfF / min(adjFits);
            bestPal = adj(idx);
            x(i,:) = x(i,:) + Mu * w * (x(bestPal(1),:) - x(i,:));
            x(i,x(i,:) > 1) = 1;
            x(i,x(i,:) < -1) = -1;
            
        end

        diff = norm(x - lastX);
        er_diffs = [er_diffs; diff];
%         if diff < Eps
%             break;
%         end

        er_xs = [er_xs; x];
    end
    % End Of CBO-ER
    
    er_meanSum(time,:) = er_means - 0.5;
    
    save('F1Partial.mat');
    
    close all;
    FName = 'F1';
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

    %FName = 'F1';
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
    
end


% ------------ === FIGURES === ------------
% %------------------------------------------------------------------------------------------------------------
save('F1FinalData.mat');

ga_meanSum = ga_meanSum - 0.5;
de_meanSum = de_meanSum - 0.5;
pso_meanSum = pso_meanSum - 0.5;
ba_meanSum = ba_meanSum - 0.5;
ws_meanSum = ws_meanSum - 0.5;
er_meanSum = er_meanSum - 0.5;

fig = figure;
FName = 'F1';
set(gcf, 'PaperPosition',[0.25 2.5 7 3]);
subplot(1,2,1);
x=-1:0.001:1;
y = - ( exp(-2*(log(2)./log(exp(1))) .* ((x-0.1)/0.8).^2) .* (sin(5*pi*x)).^6);
plot(x,y);
xlabel('\bfx');
ylabel(['\bf',FName]);
title('\bf(a)');
subplot(1,2,2);
plot(1:MaxIteration,ga_meanSum,'-*',1:MaxIteration,de_meanSum,'-^',1:MaxIteration,pso_meanSum,'-x',1:MaxIteration,ba_meanSum,'-p',1:MaxIteration,ws_meanSum,'-v',1:MaxIteration,er_meanSum,'-*');
leg = legend('GA','DE','PSO','CBO-BA','CBO-WS','CBO-ER','Location','Northeast');
set(leg,'FontSize',7);
xlabel('\bfIteration #');
ylabel(['\bf',FName]);
title('\bf(b)');
set(gca,'XScale','log');
saveas(fig,'F1.fig');
print -loos -dtiff F1.tiff;


%% Paper Figures of F1 (1 dimension)
fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 8 3]);
FName = 'F1';
subplot(1,2,1);
x = -1:0.01:1;
y = -(exp(-2.*(log(2)/log(exp(1))) .* ((x-0.1)/0.8).^2) .* (sin(5*pi.*x)).^6);
plot(x,y);
ylabel(['\bf',FName]);
xlabel('\bfx');
title('\bf(a)');
subplot(1,2,2);
plot(1:MaxIteration,ppp,'-o',1:MaxIteration,uuu,'-^',1:MaxIteration,mean(ga_meanSum),'-s',1:MaxIteration,mean(de_meanSum),'-d',1:MaxIteration,mean(pso_meanSum),'-x',1:MaxIteration,mean(ba_meanSum),'-p',1:MaxIteration,mean(ws_meanSum),'-v',1:MaxIteration,mean(er_meanSum),'-*');
title('\bf(b)');
leg = legend('DDA-C','DDA-E','GA','DE','PSO','CBO-BA','CBO-WS','CBO-ER','Location','Northeast');
set(leg,'FontSize',5);
xlabel('\bfIteration #');
set(gca,'XScale','log');
%set(gca,'YScale','log');
saveas(fig,'F1Figure.fig');
print -loos -dtiff F1Figure.tiff;


end

