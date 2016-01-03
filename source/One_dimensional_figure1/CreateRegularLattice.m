%Omid55
%Ring Regular Lattice Creator Code
%with N nodes and K neigbors (K must be even) in any node
%return sparse graph
function [ sp ] = CreateRegularLattice( N,K )

sp = sparse(N,N);
for i=1:N
    for j=1:K/2
        n1 = i + j;
        if(n1 > N)
            n1 = mod(n1 , N);
        end
        n2 = i - j;
        if(n2 <= 0)
            n2 = n2 + N;
        end
        sp(i,n1) = 1;
        sp(n1,i) = 1;
        sp(i,n2) = 1;
        sp(n2,i) = 1;
    end
end

end

