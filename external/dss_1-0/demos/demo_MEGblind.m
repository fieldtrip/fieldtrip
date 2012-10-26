% This script demonstrates the use of the DSS MATLAB package in
% separation of MEG signals using blind approach (ICA)


%   See denss for description of the denoising source separation 
%   MATLAB package

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: demo_MEGblind.m,v 1.4 2005/12/07 14:04:04 jaakkos Exp $

                                
load MEG_data.mat;                        % load the data

Xaux = X(123:end,:);                      % + 5 channels for EOG, etc.
Xaux(end,:) = -Xaux(end,:);               % make R-peaks positive
X = X(1:122,:);                           % exclude the reference channels from the data

% plot some channels and the references
% plot the results
figure
plots = 10;
for i = 1 : plots-3
  subplot(plots,1,i);
  plot(X(i,:));
end
subplot(plots,1,i+1); plot(Xaux(1,:));    % blinks
subplot(plots,1,i+2); plot(Xaux(2,:));    % saccades
subplot(plots,1,plots); plot(Xaux(5,:));  % ECG

state = dss_create_state(X)               % create the state structure 
% see dss_create_state for the definition of the different fields 
state.verbose = 3;                        % maximal information output
state.wdim = 50;                          % drop dimension to avoid overfitting
state = dss_preprocess(state.X,state);    % let's whiten for further use

state.sdim = 20;                          % extract 20 components

state = denss(state)                      % run DSS

% plot the results
figure
for i = 1 : state.sdim
  subplot(state.sdim,1,i);
  plot(state.S(i,:));
end
