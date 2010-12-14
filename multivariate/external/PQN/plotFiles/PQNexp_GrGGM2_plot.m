function [] = PQNexp_GrGGM2_plot(doBigData,nFolds)

fileName = sprintf('PQNexp_GrGGM2_%d_%d.mat',doBigData,nFolds);
load(fileName);

figure(8);
prettyPlot(lambdaVec,...
    [mean(resultL2,1);mean(resultLInf,1);mean(resultL1,1);repmat(mean(resultTik),[1 length(lambdaVec)])]',...
    {'L_{1,2}','L_{1,\infty}','L_1','Base'},...
    'Group GGM test set results','Regularization Strength (\lambda)','Average Log-Likelihood',1);
xlim([lambdaVec(1) lambdaVec(end)]);
set(gca,'Ygrid','on');