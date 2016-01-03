%Omid55
%Optimization With Opinion Formation Models
function [  ] = NodesNumber(  )

%% Clear All & Close All
clc;
close all;


%% Initiallization 
times = 50;
M = 20;
Mu = 0.9;
%Pa = 1;
MaxIteration = 10;
Eps = 10^-50;
rangeBegin = -100;
rangeEnd = 100;
optimal = rangeBegin + (rangeEnd - rangeBegin) * rand(1,M);

Nums = 15;
Means = zeros(Nums,times);
for time=1:times
    time
    for nn=1:Nums
        N = nn * 100;
        
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


     %% Calculating Values
        x = rangeBegin + (rangeEnd - rangeBegin) * rand(N,M);

        for l=1:MaxIteration
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

            diff = norm(x - lastX);
            if diff < Eps
                break;
            end
        end

        Means(nn,time) = ObjectiveFunction(mean(x),optimal) - 1;
    end
    
    close all;
    fig = figure;
    set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
    errorbar(100*(1:Nums),mean(Means'),std(Means')/sqrt(times));
    xlabel('\bfSize Of Population');
    ylabel('\bfValue');
    saveas(fig,'NodeNumber.fig');
    print -loos -dtiff NodeNumber.tiff;

    save('NumsDataPartial.mat');

end


Means = Means';

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
errorbar(100*(1:Nums),mean(Means),std(Means)/sqrt(times));
xlabel('\bfSize Of Population');
ylabel('\bfMean Cost in Last Iteration');
saveas(fig,'NodeNumber.fig');
print -loos -dtiff NodeNumber.tiff;

save('NumsDataFinal.mat');


end

