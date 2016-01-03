%Omid55
function [  ] = MuTuning(  )


%% Clear All & Close All
clc;
close all;
addpath('D:\Omid\matlab_bgl-4.0.1_2\matlab_bgl');
addpath('GA');
addpath('DE');
addpath('Dual Avg');


%% Parameters Initiallization
times = 50;
N = 100;
M = 20;
MaxIteration = 500;
Eps = 10^-50;
rangeBegin = -5;
rangeEnd = 5;


%% Runing
% Network creation
len = 10;
MuFitMeans = zeros(times,len);
for time=1:times
    time
    for Mu=0.1:0.1:1
        Mu
        averageDegree = 12;
        net = BarabasiGraphCreator(N,averageDegree/2,averageDegree/2);
        optimal = rangeBegin + (rangeEnd - rangeBegin) * rand(1,M);
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

            diff = norm(x - lastX);
            if diff < Eps
                break;
            end
        end

        MuFitMeans(time,floor(Mu/0.1)) = ObjectiveFunction(mean(x),optimal) - 1;
    end
    
    close all;
    fig = figure;
    set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
    errorbar(0.1:0.1:1,mean(MuFitMeans),std(MuFitMeans)/sqrt(time));
    xlabel('\bf\mu');
    ylabel('\bfMean Cost in Last Iteration');
    saveas(fig,'MuTunning.fig');
    print -loos -dtiff MuTunning.tiff;

    save('MuTunningDataPartial.mat');
    
end

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
errorbar(0.1:0.1:1,mean(MuFitMeans),std(MuFitMeans)/sqrt(times));
xlabel('\bf\mu');
ylabel('\bfF2');
%set(gca,'XScale','log');
%set(gca,'YScale','log');
saveas(fig,'MuTunning.fig');
print -loos -dtiff MuTunning.tiff;

save('MuTunningDataFinal.mat');

end

