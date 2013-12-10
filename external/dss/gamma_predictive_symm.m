function [params, gamma] = gamma_predictive_symm(params, state)
% Predictive controller learning rate modification rule
% for symmetric DSS
%   [gamma, params] = gamma_predictive(params, state)
%     params  Gamma function parameters
%     state   DSS algorithm state
%     gamma   Result gamma coefficient

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
    params.name = 'Predictive controller (symm)';
    params.adaptive=1;
    params.approach = {'symm'};
    return;
end

if state.iteration <= 2
  params.gamma = ones(1, size(state.W, 2));
  if state.iteration == 2
    params.deltaW = state.W_old-state.W;
  end
else
  deltaW_old = params.deltaW;
  params.deltaW = state.W_old-state.W;

  params.gamma = params.gamma + ...
      sum(deltaW_old.*params.deltaW) ./ sum(deltaW_old.^2);

  if any(params.gamma < 0.5)
    params.gamma = ones(1, size(state.W, 2)) * 0.5;
  end
  
end

gamma = params.gamma;
