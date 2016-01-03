%Omid55
function [ child ] = doCrossover( c1,c2,index1,index2 )

if exist('index2','var') == 0
    % it means 2 part crossover
    child = [c1(1:index1-1) c2(index1:end)];
    
else
    % it means 3 part crossover
    child = [c1(1:index1-1) c2(index1:index2-1) c1(index2:end)];
    
end

end

