%Omid545
function [  ] = PaTunning(  )

clc;
close all;

N = 1000;
M = 5;
Mu = 0.9;
MaxIteration = 500;
Eps = 10^-50;
rangeBegin = -100;
rangeEnd = 100;
optimal = rangeBegin + (rangeEnd - rangeBegin) * rand(1,M);

% % REAL NETWORK
% network_file = 'USairport500-net.txt';
% a = importdata(network_file);
% MM = max(max(a(:,1)),max(a(:,2)));
% sp = sparse(a(:,1),a(:,2),1,MM,MM);
% net = CreateMap(sp);
% N = size(sp,1);
% % REAL NETWORK
net = BarabasiGraphCreator(N,5,5);

result = repmat(optimal + zeros(1,M),N,1);        %  << CHECK HERE >>  FOR RESULT OF OBJECTIVE FUNCTION
x = rangeBegin + (rangeEnd - rangeBegin) * rand(N,M);


% P = randperm(size(x,2));
% save('P_rand_Perm.mat','P');

%% My Algorithm's Simulation Body (Consensus Based Optimization)(CBO)
disp('CBO');
times = 10;
firstX = x;
cost = zeros(11,times);
cnt = 0;
for time=1:times
    time
    
    for Pa=0:0.1:1
        cnt = cnt + 1;
        x = firstX;
        for l=1:MaxIteration

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
            if diff < Eps
                break;
            end
        end

        cost(cnt,time) = cost(cnt,time) + norm(x-result)*100 + l;
    end
    
    cnt = 0;
end

cost = cost';

f = figure;
set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
errorbar(0:0.1:1,mean(cost),std(cost));
xlabel('\bfPa');
ylabel('\bfCost');
saveas(f,'PaTunning.fig');
print -loos -dtiff PaTunning.tiff;

save('PaTunningData.mat');


end


