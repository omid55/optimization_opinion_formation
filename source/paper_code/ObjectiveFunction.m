%Omid55
function [ fitness ] = ObjectiveFunction( x,optimal )

%Eps = 10^-4;

% %F1's paper
% % % x = -1:0.01:1;
% % % y = (exp(2.*(log(2)/log(exp(1))) .* ((x-0.1)/0.8).^2) .* (sin(5*pi.*x)).^-6);
% % % plot(x,y);
% x = x - optimal;
% fitness = 1 / (exp(-2.*(log(2)/log(exp(1))) .* ((x-0.1)/0.8).^2) .* (sin(5*pi*x)).^6);
% %fitness = sum(1 ./ (exp(-2.*(log(2)/log(exp(1))) .* ((x-0.1)/0.8).^2) .* (sin(5*pi*x)).^6));

% fitness = -1 * ( x^4 + sin((pi*x)-1) - cos((pi*x^2)-2) );

%fitness = (2.56 - ((x-0.3) ^ 2)) / 2.56;


%% CEC

% %good
% % F1
% x = x - optimal;
% D = length(x);
% if size(x,1) ~= 1
%     x = x';
% end
%     % fitness = 0;
%     % for i=1:D
%     %     fitness = fitness + ((x(i)^2) * (10^(6*(i-1)/(D-1))));
%     % end
% fitness = sum((x.^2) .* (10.^(6.*((1:D)-1)./(D-1))));
% fitness = fitness + 1;

% good
% F2
x = x - optimal;
fitness = 0;
for i=1:length(x)
    fitness = fitness + x(i)^2-10*cos(2*pi*x(i))+10;
end
fitness = fitness + 1;

% %good
% % F3
% x = x - optimal;
% D = length(x);
% fitness = -20 * exp(-0.2*sqrt(sum(x.^2)/D)) - exp(sum(cos(2*pi*x))/D) + 20 + exp(1);
% fitness = fitness + 1;

% %bad
% % F18
% x = x - optimal;
% D = length(x);
% fitness = 0;
% load('P_rand_Perm.mat');
% m = 5;
% for k=1:floor(D/m)
%     z = x(P((k-1)*m+1:k*m));
%     for i=1:size(z,2)-1
%         fitness = fitness + (z(i)^2-z(i+1))^2 + (z(i)-1)^2;
%     end
% end
% fitness = fitness + 1;

% % F19
% x = x - optimal;
% D = length(x);
% fitness = 0;
% for i=1:D
%     s = 0;
%     for j=1:i
%         s = s + x(i);
%     end
%     fitness = fitness + s^2;
% end
% fitness = fitness + 1;

% %bad
% %F20
% x = x - optimal;
% D = length(x);
% fitness = 0;
% for i=1:D-1
%     fitness = fitness + 100 * (x(i)^2-x(i+1))^2 + (x(i)-1)^2;
% end
% %fitness = fitness + 1;

% if fitness < Eps
%     fitness = Eps;
% end


end

