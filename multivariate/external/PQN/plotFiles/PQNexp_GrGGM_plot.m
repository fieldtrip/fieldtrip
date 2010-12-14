function [] = PQNexp_GrGGM_plot(doBigData,lambda,corrections)

fileName = sprintf('PQNexp_GrGGM_%d_%f_%d.mat',doBigData,lambda,corrections);
load(fileName);

figure(4);
clf; hold on
xData = 1:200;
yData{2} = -runningMin(-fPG);
yData{1} = -runningMin(-fPQN);
yData{3} = -runningMin(-fSPG);
legendStr = {'PQN','PG','SPG'};
prettyPlot(xData,yData,legendStr,'Group GGM resuls','Function Evaluations','Objective Value');

figure(5);
clf; hold on
xData = 1:1000;
yData{2} = -runningMin(-fPG);
yData{1} = -runningMin(-fPQN);
yData{3} = -runningMin(-fSPG);
legendStr = {'PQN','PG','SPG'};
prettyPlot(xData,yData,legendStr,'Group GGM results','Function Evaluations','Objective Value');