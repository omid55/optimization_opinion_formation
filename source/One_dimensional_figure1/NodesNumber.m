%Omid55
%Optimization With Opinion Formation Models
function [  ] = NodesNumber(  )


%% Clear All & Close All
clc;
close all;


%% Initiallization 
times = 10;
M = 10;
Mu = 0.9;
Pa = 0.5;
MaxIteration = 300;
Eps = 10^-50;
rangeBegin = -32;
rangeEnd = 32;
optimal = rangeBegin + (rangeEnd - rangeBegin) * rand(1,M);

Nums = 20;
myGoodnessBests = zeros(Nums,times);
%myGoodnesses = zeros(10,1);
for time=1:times
    for nn=1:Nums
        N = nn * 100;
        result = repmat(optimal + zeros(1,M),N,1);        %  << CHECK HERE >>  FOR RESULT OF OBJECTIVE FUNCTION

        %% Network

        net = BarabasiGraphCreator(N,5,5);

%         averageDegree = 30;
%         erdosProbability = averageDegree/N;
%         sp = erdos_reyni(N,erdosProbability);
%         net = CreateMap(sp);

        % sp = WattsStrogatzCreator(N,4,0.1);
        % net = CreateMap(sp);

        % sp = CreateRegularLattice(N,4);
        % net = CreateMap(sp);

        % isClique = 1;


        % P = randperm(size(x,2));
        % save('P_rand_Perm.mat','P');


        %% Calculating Values
        x = rangeBegin + (rangeEnd - rangeBegin) * rand(N,M);

        diffs = [];
        bests = [];
        means = [];
        myGoodness = [];
        myGoodnessBest = [];
        xs = [];
        for l=1:MaxIteration
            if mod(l,10)==0
                l
            end

            % just for ploting the graph
            fits = zeros(size(x,1),1);
            for i=1:size(x,1)
                fits(i) = ObjectiveFunction(x(i,:),optimal);
            end
            means = [means; mean(fits)];
            bests = [bests; min(fits)];
            bestIdx = find(fits == min(fits));
            myGoodnessBest = [myGoodnessBest; norm(x(bestIdx(1),:) - result(1,:))];
            myGoodness = [myGoodness; norm(x - result)];
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
            diffs = [diffs; diff];
            if diff < Eps
                break;
            end

            xs = [xs; x];
        end

        myGoodnessBests(nn,time) = myGoodnessBest(end);
        %myGoodnesses(nn) = myGoodnesses(nn) + myGoodness(end);
    end
end

save('NumsData.mat');

myGoodnessBests = myGoodnessBests';

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
errorbar(100*(1:Nums),mean(myGoodnessBests),std(myGoodnessBests));
xlabel('\bfSize Of Population');
ylabel('\bfValue');
saveas(fig,'NodeNumber.fig');
print -loos -dtiff NodeNumber.tiff;

% pause(30);
% exit;


end

