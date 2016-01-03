%Omid55
%NewOrder-Based Crossovers for the Graph Coloring Problem
%Christine L. Mumford's Method

function [ gradient ] = GradientOfCalculateFitness_real( x )

x = CalculateValues(x);
Ln2 = (log(2)/log(exp(1)));
gradient = -1 / (1+exp(-2*(log(2)/log(exp(1))) * ((x-0.1)/0.8).^2) * (sin(5*pi*x))^6);   %5*Ln2*((x-0.1)/0.8) * exp(2*Ln2*((x-0.1)/0.8)^2) * (sin(5*pi*x))^-6 + (-30*pi)*cos(5*pi*x)*((sin(5*pi*x))^-7) * exp(2*Ln2*((x-0.1)/0.8)^2);

end

