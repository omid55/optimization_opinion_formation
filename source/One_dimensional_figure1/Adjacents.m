function [ nodes ] = Adjacents( sp,nodeIndex )

nodes = [OuterNodes(sp,nodeIndex) InnerNodes(sp,nodeIndex)];
nodes = unique(nodes);

end

