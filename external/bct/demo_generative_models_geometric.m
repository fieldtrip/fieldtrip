% Demonstration of generative model functions. See GENERATIVE_MODEL and
% EVALUATE_GENERATIVE_MODEL for further details and interpretation.

clear
close all
clc

data = load('demo_generative_models_data');
A     = data.A;
Aseed = data.Aseed;
D     = data.D;

% get cardinality of network
n = length(A);

% set model type
modeltype = 'sptl';

% set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'},{'powerlaw'}];

% choose some model parameters
nparams = 100;
params = unifrnd(-10,0,nparams,1);

% generate synthetic networks and energy for the neighbors model;
[B,E,K] = evaluate_generative_model(Aseed,A,D,modeltype,modelvar,params);
X = [E,K];

% show scatterplot of parameter values versus energy and KS statistics
names = [...
    {'energy'},...
    {'degree'},...
    {'clustering'},...
    {'betweenness'},...
    {'edge length'}];

f = figure(...
        'units','inches',...
        'position',[2,2,4,4]);
for i = 1:size(X,2)
    subplot(3,2,i);
    scatter(params,X(:,i),100,X(:,i),'filled');
    set(gca,...
        'ylim',[0,1],...
        'clim',[0,1]);
    colormap(jet);
    xlabel('geometric parameter, \eta');
    ylabel(names{i});
end
