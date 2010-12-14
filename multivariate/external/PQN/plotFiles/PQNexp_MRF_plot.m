function [] = PQNexp_MRF_plot(doBigData,nInstances,lambda,corrections,infer)

fileName = sprintf('PQNexp_MRF_%d_%d_%f_%d_%s.mat',doBigData,nInstances,lambda,corrections,infer);
load(fileName);

figure(6);
clf; hold on
xData = 1:100;
yData{3} = runningMin(fGraft);
yData{2} = runningMin(fSPG);
yData{1} = runningMin(fPQN2);
yData{4} = runningMin(fPQN);
legendStr = {'PQN','SPG','Graft','PQN-SOC'};
prettyPlot(xData,yData,legendStr,'MRF results','Function Evaluations','Objective Value');

figure(7);
clf; hold on
xData = 1:1000;
yData{3} = runningMin(fGraft);
yData{2} = runningMin(fSPG);
yData{1} = runningMin(fPQN2);
yData{4} = runningMin(fPQN);
legendStr = {'PQN','SPG','Graft','PQN-SOC'};
prettyPlot(xData,yData,legendStr,'MRF results','Function Evaluations','Objective Value');
