%Omid55
% Get n distinct items from range 
function [ result ] = GetDistinctItems( list,n )

result = [];
for i=1:n
    index = randi(length(list),1,1);
    result = [result; list(index)];
    list(index) = [];
end

end

