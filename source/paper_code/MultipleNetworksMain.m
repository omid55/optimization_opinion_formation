%Omid55
function [  ] = MultipleNetworksMain(  )

clc;
close all;

% addpath('D:\Omid\matlab_bgl-4.0.1_2\matlab_bgl');

times = 10;
N = 1000;
M = 10;
Mu = 0.6;
Pa = 1;
MaxIteration = 400;
Eps = 10^-50;
rangeBegin = -5;
rangeEnd = 5;
optimal = rangeBegin + (rangeEnd - rangeBegin) * rand(1,M);
result = repmat(optimal + zeros(1,M),N,1);        %  << CHECK HERE >>  FOR RESULT OF OBJECTIVE FUNCTION

x = rangeBegin + (rangeEnd - rangeBegin) * rand(N,M);
firstX = x;

sumsMean = cell(5,1);
sumsBest = cell(5,1);
for time=1:times
    time
    
    myMeansCell = cell(5,1);
    myBestsCell = cell(5,1);
    for netw=1:5
        x = firstX;
        netw

        switch netw
            case 1
                net = BarabasiGraphCreator(N,4,4);

            case 2
                erdosProbability = 0.1;
                sp = erdos_reyni(N,erdosProbability);
                net = CreateMap(sp);

            case 3
                sp = WattsStrogatzCreator(N,4,0.1);
                net = CreateMap(sp);

            case 4
                sp = CreateRegularLattice(N,4);
                net = CreateMap(sp);

            case 5
                %Clique
        end
        
%         P = randperm(size(x,2));
%         save('P_rand_Perm.mat','P');

     %% My Algorithm's Simulation Body (Consensus Based Optimization)(CBO)
        %disp('CBO');

        myGoodness = [];
        myGoodnessBest = [];
        myMeans = [];
        myBests =[];
        for l=1:MaxIteration
            
            % just for ploting the graph
            fits = zeros(size(x,1),1);
            for i=1:size(x,1)
                fits(i) = ObjectiveFunction(x(i,:),optimal);
            end
            myMeans = [myMeans; mean(fits)];
            myBests = [myBests; min(fits)];
            bestIdx = find(fits == min(fits));
            myGoodnessBest = [myGoodnessBest; norm(x(bestIdx(1),:) - result(1,:))];
            myGoodness = [myGoodness; norm(x - result)];
            % just for plotting the graph

            lastX = x;

            for i=1:N  % for each agent
                if net == 5    %Clique
                    adj = 1:N;
                    adj(i) = [];
                else
                    adj = net(num2str(i));
                end
                
                if length(adj) == 0
                    break;
                end

                dummy = 0.2;      % << CHECK HERE >>   FOR VALUE OF DUMMY VARIABLE
                for feat=1:M
                    fits = zeros(length(adj),1);
                    for j=1:length(adj)
                        mask = dummy * ones(1,M);   % dummy values
                        mask(feat) = x(adj(j),feat);
                        fits(j) = ObjectiveFunction(mask,optimal);
                    end
                    bestIdx = find(fits == min(fits));
                    bestAdj = adj(bestIdx(1));
                    selfMask = dummy * ones(1,M);
                    selfMask(feat) = x(i,feat);
                    selfMaskedFit = ObjectiveFunction(selfMask,optimal);
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
            if diff < Eps
                break;
            end
        end
        
        myMeansCell{netw} = myMeans;
        myBestsCell{netw} = myBests;
    end
    
    for i=1:size(sumsMean)
        sumsMean{i} = AddUnmatchSizeVectors(sumsMean{i},myMeansCell{i})';
        sumsMean{i} = AddUnmatchSizeVectors(sumsMean{i},myMeansCell{i})';
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:size(sumsMean)
        myMeansCell{i} = sumsMean{i}/time;
    end

    for i=1:length(myMeansCell)
        myMeansCell{i} = myMeansCell{i} - 1;
        vec = ones(MaxIteration,1) * myMeansCell{i}(end);
        vec(1:length(myMeansCell{i})) = myMeansCell{i};
        myMeansCell{i} = vec;
    end

    close all;
    fig = figure;
    set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
    plot(1:length(myMeansCell{1}),myMeansCell{1},'-v',1:length(myMeansCell{2}),myMeansCell{2},'-*',1:length(myMeansCell{3}),myMeansCell{3},'-^',1:length(myMeansCell{4}),myMeansCell{4},'-d',1:length(myMeansCell{5}),myMeansCell{5},'-p');
    leg = legend('BA','ER','WS','Regular','Clique','Location','Northeast');
    xlabel('\bfIteration #');
    ylabel('\bfF2');
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    %title('\bfDifference From Optimal Result');
    saveas(fig,'Nets(Mean).fig');
    print -loos -dtiff Nets(Mean).tiff;
    fig = figure;
    set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
    plot(1:length(myMeansCell{1}),myMeansCell{1},'-v',1:length(myMeansCell{2}),myMeansCell{2},'-*',1:length(myMeansCell{3}),myMeansCell{3},'-^',1:length(myMeansCell{4}),myMeansCell{4},'-d',1:length(myMeansCell{5}),myMeansCell{5},'-p');
    leg = legend(leg,'BA','ER','WS','Regular','Clique','Location','Northeast');
    set(leg,'FontSize',7);
    xlabel('\bfIteration #');
    ylabel('\bfF2');
    set(gca,'XScale','log');
    %title('\bfDifference From Optimal Result');
    saveas(fig,'LogNets(Mean).fig');
    print -loos -dtiff LogNets(Mean).tiff;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

for i=1:size(sumsMean)
    myMeansCell{i} = sumsMean{i}/times;
end

for i=1:length(myMeansCell)
    myMeansCell{i} = myMeansCell{i} - 1;
    vec = ones(MaxIteration,1) * myMeansCell{i}(end);
    vec(1:length(myMeansCell{i})) = myMeansCell{i};
    myMeansCell{i} = vec;
end


fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
plot(1:length(myMeansCell{1}),myMeansCell{1},'-v',1:length(myMeansCell{2}),myMeansCell{2},'-*',1:length(myMeansCell{3}),myMeansCell{3},'-^',1:length(myMeansCell{4}),myMeansCell{4},'-d',1:length(myMeansCell{5}),myMeansCell{5},'-p');
leg = legend('BA','ER','WS','Regular','Clique','Location','Northeast');
xlabel('\bfIteration #');
ylabel('\bfF2');
set(gca,'XScale','log');
set(gca,'YScale','log');
%title('\bfDifference From Optimal Result');
saveas(fig,'Nets(Mean).fig');
print -loos -dtiff Nets(Mean).tiff;

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
plot(1:length(myMeansCell{1}),myMeansCell{1},'-v',1:length(myMeansCell{2}),myMeansCell{2},'-*',1:length(myMeansCell{3}),myMeansCell{3},'-^',1:length(myMeansCell{4}),myMeansCell{4},'-d',1:length(myMeansCell{5}),myMeansCell{5},'-p');
leg = legend(leg,'BA','ER','WS','Regular','Clique','Location','Northeast');
set(leg,'FontSize',7);
xlabel('\bfIteration #');
ylabel('\bfF2');
set(gca,'XScale','log');
%title('\bfDifference From Optimal Result');
saveas(fig,'LogNets(Mean).fig');
print -loos -dtiff LogNets(Mean).tiff;

% fig = figure;
% set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
% plot(1:length(myBestsCell{1}),myBestsCell{1},'-v',1:length(myBestsCell{2}),myBestsCell{2},'-*',1:length(myBestsCell{3}),myBestsCell{3},'-^',1:length(myBestsCell{4}),myBestsCell{4},'-d',1:length(myBestsCell{5}),myBestsCell{5},'-p');
% leg = legend(leg,'BA','ER','WS','Regular','Clique','Location','Northeast');
% xlabel('\bfIteration #');
% ylabel('\bfF2');
% %title('\bfDifference From Optimal Result');
% saveas(fig,'Nets(Best).fig');
% print -loos -dtiff Nets(Best).tiff;
% 
% fig = figure;
% set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
% plot(1:length(myBestsCell{1}),myBestsCell{1},'-v',1:length(myBestsCell{2}),myBestsCell{2},'-*',1:length(myBestsCell{3}),myBestsCell{3},'-^',1:length(myBestsCell{4}),myBestsCell{4},'-d',1:length(myBestsCell{5}),myBestsCell{5},'-p');
% leg = legend(leg,'BA','ER','WS','Regular','Clique','Location','Northeast');
% set(leg,'FontSize',7);
% xlabel('\bfIteration #');
% ylabel('\bfF2');
% set(gca,'XScale','log');
% %title('\bfDifference From Optimal Result');
% saveas(fig,'LogNets(Best).fig');
% print -loos -dtiff LogNets(Best).tiff;


save('NetsData.mat');


end

