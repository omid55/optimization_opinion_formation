index = 1;
xss = cell(size(xs,1)/N,1);
for i=1:size(xs,1)/N
    xss{i} = xs(index:index+N-1,:);
    index = index+N;
end

f = figure;
set(gcf, 'PaperPosition',[0.25 2.5 7 3.5]);
subplot(1,2,1);
for i=1:N
    ags = zeros(size(xss,1),M);
    for j=1:size(xss,1)
        ab = xss{j};
        ab = ab(i,:);
        ags(j,:) = ab;
    end
    plot(ags);
    hold on;
end
%title('\bfOpinion Changing');
title('\bf(a)');
xlabel('\bfIteration #');
ylabel('\bfOpinion');

subplot(1,2,2),plot(diffs);
%title('\bfGrow of Opinions');
title('\bf(b)');
xlabel('\bfIteration #');
ylabel('\bfDifference');
saveas(f,'OpinionsChanges.fig');
print -loos -dtiff OpinionsChanges.tiff;