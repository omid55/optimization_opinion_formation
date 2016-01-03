%Omid55
function [ sp ] = WattsStrogatzCreator( N,K,P )

probability = 1 - P;
sp = CreateRegularLattice(N,K);
for i=1:N
    adj = Adjacents(sp,i);
    for j=1:length(adj)
        if(rand(1)>probability)
            available = setdiff(1:N,[adj i]); 
            index = 1 + floor(rand(1)*length(available));
            sp(i,adj(j)) = 0;
            sp(adj(j),i) = 0;
            sp(i,available(index)) = 1;
            sp(available(index),i) = 1;
        end
    end
end

end

