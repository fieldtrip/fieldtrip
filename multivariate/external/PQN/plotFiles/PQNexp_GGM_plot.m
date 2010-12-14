function [] = PQNexp_GGM_plot(doBigData,lambda,corrections)

fileName = sprintf('PQNexp_GGM_%d_%f_%d.mat',doBigData,lambda,corrections);
load(fileName);

figure(2);
clf; hold on
xData = 1:200;
yData{3} = -runningMin(-fBCD);
yData{2} = -runningMin(-fPG);
yData{1} = -runningMin(-fPQN);
yData{4} = -runningMin(-fSPG);
legendStr = {'PQN','PG','BCD','SPG'};
prettyPlot(xData,yData,legendStr,'GGM results','Function Evaluations','Objective Value');

figure(3);
clf; hold on
xData = 1:1000;
prettyPlot(xData,yData,legendStr,'GGM results','Function Evaluations','Objective Value');