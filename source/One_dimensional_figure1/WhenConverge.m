%Omid55
% this function tell us when array of output of the algorithm converged and
% in which index.
function [ index ] = WhenConverge( array )

index = length(array);   % this means if not converged the function will return maxIteration
Eps = 10^-3;
for i=1:length(array)-1
    if(abs(array(i)-array(i+1)) < Eps)
        index = i;
        break;
    end
end


end

