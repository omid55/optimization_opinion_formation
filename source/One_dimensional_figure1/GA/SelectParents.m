%Omid55
function [ parents ] = SelectParents ( fitnesses )

disp('SelectParents ... ');
% Stochastic Universal Sampling
mu = length(fitnesses);
parents = [];
r = (1/mu) * rand(1,1);

for i=1:mu
    
    s = 0;
    j = 0;
    
    while s <= r
        j = j + 1;
        s = s + fitnesses(j);
    end
    
    parents = [parents j];
    r = r + (1/mu);
    
end

end