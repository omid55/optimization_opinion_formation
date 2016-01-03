%Omid55
function [ nodes ] = InnerNodes( sp,nodeIndex )

nodes = [];
for i=1:size(sp,1)
    if(sp(i,nodeIndex)==1)
        nodes = [nodes i];
    end
end

end

