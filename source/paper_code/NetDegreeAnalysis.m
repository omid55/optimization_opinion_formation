%Network Degree Analysis
degs = zeros(N,1);
for i=1:N
	degs(i) = length(net(num2str(i)));
end
hist(degs);
mean(degs)
min(degs)
max(degs)
