%Omid55
function [ nodes ] = OuterNodes( sp,nodeIndex )

nodes = [];
for i=1:size(sp,1)
    if(sp(nodeIndex,i)==1)
        nodes = [nodes i];
    end
end

end

