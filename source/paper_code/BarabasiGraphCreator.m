%Omid55
function [ net ] = BarabasiGraphCreator( N,m,m0 )

[edges nodes] = CliqueCreator(m0);
net = containers.Map;
for i=1:size(edges,1)
    node = num2str(edges(i,1));
    if net.isKey(node) == 0
        net(node) = edges(i,2);
    else
        net(node) = [net(node) edges(i,2)];
    end
end
probablity = repmat(nodes,1,m0);

n = m0;
while n < N
    n = n + 1;
    selected = [];
    while length(selected) < m
        if(length(probablity) == 0)
            selected = 1;
        else
            prob = probablity;
            for i=1:length(selected)
                inds = find(prob == selected(i));
                prob(inds) = [];
            end
            index = 1 + floor(rand(1)*length(prob));
            selected = [selected; prob(index)];
            selected = unique(selected);   %
        end
    end
    selected = selected';
    nStr = num2str(n);
    if net.isKey(nStr) == 0
        net(nStr) = selected;
    else
        net(nStr) = [net(nStr) selected];
    end
    for i=1:length(selected)
        node = num2str(selected(i));
        if net.isKey(node) == 0
            net(node) = n;
        else
            net(node) = [net(node) n];
        end
    end
    probablity=[probablity selected ones(1,length(selected))*n];
end

end

