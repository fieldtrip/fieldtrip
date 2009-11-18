% This Matlab script demonstrates how to use the neuroML toolbox
% together with FieldTrip in order to classify EEG/MEG data based on
% single trials. The dataset is already preprocessed and freqanalysed and
% consists of a cell array representing left motor execution (Lx), left motor
% imagery (Li), right motor execution (Rx), right motor imagery (Ri), and no motor
% response (No) classes.

%% clear data and fix the random number generator (RNG)

fclose('all');
close all
clear all

format long

% fixing the RNG in order to reproduce the experiment
rand('twister',1); randn('state',1);

% create model data

factors = cell(1,2);
factors{1} = multinomial_cpd(1,[],[0.5 0.5]');
factors{2} = gaussian_cpd(2,[],1,[-0.01 0.01]',cell(2,1),[0.001 0.001]');

bn = bayesnet(factors);

% simulate data
data = bn.simulate(5000);

clf = rnb('lambda',0);
clf = clf.train(data(:,2),data(:,1));

clf2 = nb();
clf2 = clf2.train(data(:,2),data(:,1));

