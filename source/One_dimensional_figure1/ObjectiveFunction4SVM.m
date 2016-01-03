%Omid55
function [ fitness ] = ObjectiveFunction4SVM( x,dataTrain,labelsTrain,dataTest,labelsTest )

Eps = 10^-30;

%% SVM Parameter Tuning
%load('MyData.mat');
while x(4)>90
    x(4) = x(4) - 10;
end
params = strcat(['-s ' num2str(floor(mod(x(1),2))) ' -t ' num2str(floor(mod(x(2),4))) ' -r 1 -c ' num2str(abs(x(3))) ' -n ' num2str(abs(x(4))) ' -b 1 -q']);
svmModel = svmtrain(labelsTrain,dataTrain,params);
[pred,acc] = svmpredict(zeros(size(dataTest,1),1),dataTest,svmModel,'-b 1');
if size(pred,1)==0
    fitness = 1;
else
    error = length(find(pred ~= labelsTest))/length(labelsTest);
    fitness = error;   % minimizing error
end

% %% NN Parameter Tuning
% %load('MyData.mat');
% x=abs(x);
% numLayers = ceil(x(1) * 3);
% numNeurons = ceil(x(2) * 5);
% LR = x(3);
% if numLayers > 3
%     numLayers = 3;
% end
% if numNeurons > 5
%     numNeurons = 5;
% else
%     if numNeurons == 0
%         numNeurons = 1;
%     end
% end
% if LR > 1
%     LR = 1;
% end
% layers=ones(1,numLayers)*numNeurons;  
% net=newff(dataTrain',labelsTrain',layers,{},'trainlm');
% net.trainParam.lr=LR;
% net.trainParam.mc=0.9;
% net.trainParam.show=25;
% net.trainParam.epochs=1000;
% net.trainParam.goal=0;
% net.trainParam.max_fail=6;
% net.trainParam.min_grad=1e-010;
% net.trainParam.time=Inf;
% net.trainParam.mu=0.001;
% net.trainParam.mu_dec=0.1;
% net.trainParam.mu_inc=10;
% net.trainParam.mu_max=10000000000;
% %net.trainParam.mem_reduc=1;
% net.trainParam.showCommandLine=false;
% net.trainParam.showWindow=false;
% net=train(net,dataTrain',labelsTrain');
% res=sim(net,dataTest');
% fitness = mse(res'-labelsTest);

if fitness == 0
    fitness = Eps;
end

end

