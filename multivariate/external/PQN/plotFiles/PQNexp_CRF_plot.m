function [] = PQNexp_CRF_plot(doBigData,lambda,nInstances,nVars)

fileName = sprintf('PQNexp_CRF_%d_%d_%d_%d.mat',doBigData,lambda,nInstances,nVars);
load(fileName);

figure(1);
clf;clear yData
maxPlotIter = 500;
xData = 1:maxPlotIter;
yData{1} = runningMin(fPQN);
yData{2} = runningMin(fSPG);
yData{3} = runningMin(fOW);
yData{4} = runningMin(fPr);
yData{5} = runningMin(fPQNB);
legendStr = {'PQN','SPG','OWL-QN','L-BFGS-B','PQN-B'};
prettyPlot(xData,yData,legendStr,'CRF results','Function Evaluations','Objective Value');