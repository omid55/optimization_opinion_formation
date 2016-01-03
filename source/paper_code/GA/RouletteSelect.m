%Omid55
function [ parent ] = RouletteSelect( fitnesses )

r = rand(1,1) * sum(fitnesses);
s = 0;
j = 0;

while s <= r
    j = j + 1;
    s = s + fitnesses(j);
end

parent = j;


end

