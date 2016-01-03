%Omid55
function [  ] = MultipleNetworksMain(  )

clc;
close all;

addpath('D:\Omid\matlab_bgl-4.0.1_2\matlab_bgl');

times = 1;
N = 5000;
M = 10;
Mu = 0.9;
Pa = 1;
MaxIteration = 400;
Eps = 10^-50;
rangeBegin = -100;
rangeEnd = 100;
optimal = rangeBegin + (rangeEnd - rangeBegin) * rand(1,M);
result = repmat(optimal + zeros(1,M),N,1);        %  << CHECK HERE >>  FOR RESULT OF OBJECTIVE FUNCTION

x = rangeBegin + (rangeEnd - rangeBegin) * rand(N,M);
firstX = x;

sums = cell(5,1);
for time=1:times
    time
    
    myGoodnesses = cell(5,1);
    for netw=1:5
        x = firstX;
        netw

        switch netw
            case 1
                net = BarabasiGraphCreator(N,5,5);

            case 2
                averageDegree = 30;
                erdosProbability = averageDegree/N;
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
        
        P = randperm(size(x,2));
        save('P_rand_Perm.mat','P');

        %% My Algorithm's Simulation Body (Consensus Based Optimization)(CBO)
        disp('CBO');

        myGoodness = [];
        myGoodnessBest = [];
        for l=1:MaxIteration
            if mod(l,10)==0
                l
            end

            % just for ploting the graph
            fits = zeros(size(x,1),1);
            for i=1:size(x,1)
                fits(i) = ObjectiveFunction(x(i,:),optimal);
            end
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
        
        myGoodnesses{netw} = myGoodness;
    end
    
    for i=1:size(myGoodnesses)
        sums{i} = AddUnmatchSizeVectors(sums{i},myGoodnesses{i})';
    end
    
end

for i=1:size(myGoodnesses)
    myGoodnesses{i} = sums{i}/times;
end

fig = figure;
leg = plot(1:length(myGoodnesses{1}),myGoodnesses{1},'-v',1:length(myGoodnesses{2}),myGoodnesses{2},'-*',1:length(myGoodnesses{3}),myGoodnesses{3},'-^',1:length(myGoodnesses{4}),myGoodnesses{4},'-d',1:length(myGoodnesses{5}),myGoodnesses{5},'-p');
legend(leg,'BA','ER','WS','Regular','Clique','Location','best');
xlabel('\bfIteration #');
ylabel('\bfDifference Norm');
%title('\bfDifference From Optimal Result');
saveas(fig,'Nets.fig');
print -loos -dtiff Nets.tiff;

fig = figure;
leg = plot(1:length(myGoodnesses{1}),log(myGoodnesses{1}),'-v',1:length(myGoodnesses{2}),log(myGoodnesses{2}),'-*',1:length(myGoodnesses{3}),log(myGoodnesses{3}),'-^',1:length(myGoodnesses{4}),log(myGoodnesses{4}),'-d',1:length(myGoodnesses{5}),log(myGoodnesses{5}),'-p');
legend(leg,'BA','ER','WS','Regular','Clique','Location','best');
xlabel('\bfIteration #');
ylabel('\bfDifference Norm');
%title('\bfDifference From Optimal Result');
saveas(fig,'LogNets.fig');
print -loos -dtiff LogNets.tiff;

save('NetsData.mat');


end

