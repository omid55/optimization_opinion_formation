%Omid55
function [ fitness ] = CalculateFitness_real( x )

xv = CalculateValues(x);
fitness = 1 / (1+exp(-2*(log(2)/log(exp(1))) * ((xv-0.1)/0.8).^2) * (sin(5*pi*xv))^6);

end

