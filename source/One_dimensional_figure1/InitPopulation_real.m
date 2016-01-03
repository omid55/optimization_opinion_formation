%Omid55
function [ population ] = InitPopulation_real( x ) 

%disp('InitPopulation ... ');
PrecissionRate = 5;
pop = dec2bin(floor(abs(x) * 10^PrecissionRate));
population = zeros(size(pop,1),size(pop,2)+1);
population(:,1) = sign(x);
for i=1:size(pop,2)
    population(:,i+1) = str2num(pop(:,i));
end


end