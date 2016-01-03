%Omid55
function [ gradient ] = GradientOfObjectiveFunction( x,optimal,d )

%% Gradient of CEC functions for each dimension

% %F1
% x = x - optimal;
% D = length(x);
% gradient = 2 * x(d) * (10 ^ (6*(d-1)/(D-1)));

%F2
x = x - optimal;
gradient = 2*x(d) + 20*pi*sin(2*pi*x(d));

% %F3
% x = x - optimal;
% D = length(x);
% p1 = sum(x.^2);
% p2 = sum(cos(2*pi*x));
% gradient = (4/D)*x(d)*((p1/D)^-0.5)*exp(-0.2*sqrt(p1/D)) - (1/D)*-2*pi*sin(2*pi*x(d))*exp(p2/D);

% %F19
% x = x - optimal;
% D = length(x);
% gradient = 0;
% for i=d:D
%     s = 0;
%     for j=1:i
%         s = s + x(i);
%     end
%     gradient = gradient + s * 2;
% end


% %F1's paper
% x = x - optimal;
% Ln2 = (log(2)/log(exp(1)));
% gradient = 5*Ln2*((x-0.1)/0.8) * exp(2*Ln2*((x-0.1)/0.8)^2)*(sin(5*pi*x))^-6 + (-30*pi)*cos(5*pi*x)*((sin(5*pi*x))^-7) * exp(2*Ln2*((x-0.1)/0.8)^2);
% %gradient = sum(5*Ln2.*((x-0.1)./0.8).*exp(2*Ln2*((x-0.1)./0.8).^2).*(sin(5*pi*x)).^-6 + (-30*pi)*cos(5*pi*x).*((sin(5*pi*x)).^-7) .* exp(2*Ln2.*((x-0.1)./0.8).^2));


end

