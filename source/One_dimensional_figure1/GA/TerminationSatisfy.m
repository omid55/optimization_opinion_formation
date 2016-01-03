%Omid55
function [ terminate ] = TerminationSatisfy( fitnesses,Eps )

Best = min(fitnesses)
Mean = mean(fitnesses) 

leng = length(find(fitnesses <= Mean));
if leng >= 1 * length(fitnesses) || abs(Best - Mean)<Eps     %if leng >= 0.9 * length(fitnesses) || abs(Best - Mean)<Eps
    terminate = 1;
else
    terminate = 0;
end

end

