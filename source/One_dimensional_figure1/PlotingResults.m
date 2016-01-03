%Omid55

clear;
close all;
clc;

% %% Simple Function
% load('BA-simple.mat');
% ba_diffs = diffs;
% load('WS-simple.mat');
% ws_diffs = diffs;
% fig = figure;
% xx = -1:0.01:1;
% yy = - ( exp(-2*(log(2)/log(exp(1))) .* ((xx-0.1)/0.8).^2) .* (sin(5*pi*xx)).^6);
% subplot(1,2,1),plot(xx,yy);
% xlabel('\bfInput');
% ylabel('\bfOutput');
% title('\bf(a)');
% subplot(1,2,2),set(gcf, 'PaperPosition',[0.25 2.5 8 3.5]);
% plot(1:size(ba_diffs,1),ba_diffs,'-v',1:size(ws_diffs,1),ws_diffs,'-*');
% legend('BA','WS','Location','Northeast');
% xlabel('\bfIteration #');
% ylabel('\bfDifference');
% title('\bf(b)');
% saveas(fig,'SimpleFunction.fig');
% print -loos -dtiff SimpleFunction.tiff;
% %%

load('BA-TheWholeData.mat');
ba_ga_goodness = ga_goodness;
ba_ga_bests = ga_bests;
ba_ga_means = ga_means;
ba_ga_generation = ga_generation;
ba_ga_goodnessBest = ga_goodnessBest;
ba_de_goodness = de_goodness;
ba_de_goodnessBest = de_goodnessBest;
ba_de_bests = de_bests;
ba_de_means = de_means;
ba_de_it = de_it;
ba_pso_goodness = pso_goodness;
ba_pso_goodnessBest = pso_goodnessBest;
ba_pso_bests = pso_bests;
ba_pso_means = pso_means;
ba_pso_it = pso_it;
ba_xs = xs;
ba_myGoodness = myGoodness;
ba_myGoodnessBest = myGoodnessBest;
ba_diffs = diffs;
ba_bests = bests;
ba_means = means;
ba_LastIteration = LastIteration;

load('WS-TheWholeData.mat');
ws_ga_goodness = ga_goodness;
ws_ga_goodnessBest = ga_goodnessBest;
ws_ga_bests = ga_bests;
ws_ga_means = ga_means;
ws_ga_generation = ga_generation;
ws_de_goodness = de_goodness;
ws_de_goodnessBest = de_goodnessBest;
ws_de_bests = de_bests;
ws_de_means = de_means;
ws_de_it = de_it;
ws_pso_goodness = pso_goodness;
ws_pso_goodnessBest = pso_goodnessBest;
ws_pso_bests = pso_bests;
ws_pso_means = pso_means;
ws_pso_it = pso_it;
ws_xs = xs;
ws_myGoodness = myGoodness;
ws_myGoodnessBest = myGoodnessBest;
ws_diffs = diffs;
ws_bests = bests;
ws_means = means;
ws_LastIteration = LastIteration;


%% Figures

%------------------------------------------------------------------------------------------------------------
fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 6 5]);

subplot(2,2,1),plot(1:size(ba_diffs,1),ba_diffs,'-v',1:size(ws_diffs,1),ws_diffs,'-*');
set(gca,'XScale','log');
legend('BA','WS','Location','Northeast');
xlabel('\bfIteration #');
ylabel('\bfDifference');
title('\bf(a)');

subplot(2,2,2),plot(1:ba_ga_generation,ba_ga_means,'-*',1:ba_de_it,ba_de_means,'-^',1:ba_pso_it,ba_pso_means,'-x',1:ba_LastIteration,ba_means,'-p',1:ws_LastIteration,ws_means,'-o');
set(gca,'XScale','log');
ylabel('\bfValue');
xlabel('\bfIteration #');
title('\bf(b)');
legend('GA','DE','PSO','CBO-BA','CBO-WS','Location','East');

subplot(2,2,3),plot(1:ba_ga_generation,ba_ga_goodness,'-v',1:ba_de_it,ba_de_goodness,'-^',1:ba_pso_it,ba_pso_goodness,'-x',1:ba_LastIteration,ba_myGoodness,'-*',1:ws_LastIteration,ws_myGoodness,'-o');
set(gca,'XScale','log');
xlabel('\bfIteration #');
ylabel('\bfDifference Norm');
title('\bf(c)');

subplot(2,2,4), plot(1:ba_ga_generation,ba_ga_goodnessBest,'-v',1:ba_de_it,ba_de_goodnessBest,'-^',1:ba_pso_it,ba_pso_goodnessBest,'-x',1:ba_LastIteration,ba_myGoodnessBest,'-*',1:ws_LastIteration,ws_myGoodnessBest,'-o');
set(gca,'XScale','log');
xlabel('\bfIteration #');
ylabel('\bfDifference Norm');
title('\bf(d)');

saveas(fig,'MyFigure.fig');
print -loos -dtiff MyFigure.tiff;
%------------------------------------------------------------------------------------------------------------

