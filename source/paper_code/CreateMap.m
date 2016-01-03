%Omid55
function [ mapNet ] = CreateMap( netSp )

mapNet = containers.Map;
for i=1:size(netSp,1)
    adj = Adjacents(netSp,i);
    mapNet(num2str(i)) = adj;
end

end

