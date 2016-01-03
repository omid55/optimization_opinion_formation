%Omid55
function [ xv ] = CalculateValues( xs )

xv = zeros(size(xs,1),1);
for i=1:size(xs,1)
    x = xs(i,:);
    PrecissionRate = 5;
    sgn = sign(x(1));
    x = x(2:end);
    for j=1:length(x)
        while x(j) > 1
            x(j) = x(j) - 1;
        end
        while x(j) < -1
            x(j) = x(j) + 1;
        end
    end

    xv(i) = sgn * bin2dec(num2str(round(abs(x))))/10^PrecissionRate;
end

end