%Omid55
function [ c ] = AddUnmatchSizeVectors( a,b )

if size(a,1) ~= 1
    a = a';
end
if size(b,1) ~= 1
    b = b';
end

if length(a) > length(b)
    c = a + [b  zeros(1, length(a) - length(b))];
else
    c = b + [a  zeros(1, length(b) - length(a))];
end

end


