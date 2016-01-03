%Omid55
function [ edges,nodes ] = CliqueCreator( k )

edges = [];
nodes = 1:k;
for i=1:k
    for j=1:k
        if(i~=j)
            edges = [edges;[i j]]; 
        end
    end
end

end

