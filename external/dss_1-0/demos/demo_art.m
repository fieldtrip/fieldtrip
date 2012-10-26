% This script demonstrates the use of the DSS MATLAB package in
% separation of artificially generated signals


%   See denss for description of the denoising source separation 
%   MATLAB package

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: demo_art.m,v 1.3 2005/12/07 11:24:14 jaakkos Exp $


load art_data                % load the mixtures
 

state = dss_create_state(X)  % create the state structure 

% state contains (as the default) the following fields:
% see help dss_create_state for more information

% X: [5x8192 double]      the data
% verbose: 1              amount of output (minimal)
% algorithm: 'defl'       deflation (iterative DSS extracting one component at a time) 
% preprocf: [1x1 struct]  default preprocessing function (PCA sphering)
% orthof: [1x1 struct]    orthogonalisation function (ortho_default)
% denf: [1x1 struct]      denoising function (FastICA tanh as the default)
% stopf: [1x1 struct]     stopping criterion function (default_stop)
% report: [1x1 struct]    reports that are gathered per iteration or per component 
%                         (default: none, see report_*.m functions for possibilities)
%
% alpha: 1                scaling constant for the denoising
% beta: 0                 spectral shift (by default: none)
%                         the final denoising will be s_den = alpha * (denf(s) + beta * s)
% gamma: 1                step size (by default 1) 
%                         see gamma_179 and gamma_predictive for alternatives)


state.verbose = 3;        % maximal information output
                          % other fields may be changed like this too


state = denss(state)      % run the dss algorithm

% plot the results
figure
for i = 1 : state.sdim
  subplot(state.sdim,1,i);
  plot(state.S(i,:));
end
