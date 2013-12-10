function [params, gamma] = gamma_predictive(params, state)
% Predictive controller learning rate modification rule
% for deflation DSS
%   [params, gamma] = gamma_predictive(params, state)
%     params  Gamma function parameters
%     state   DSS algorithm state
%     gamma   Result gamma coefficient

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
    params.name = 'Predictive controller';
    params.adaptive=1;
    params.approach = {'defl'};
    return;
end

if state.iteration <= 2
  params.gamma = 1;
  if state.iteration == 2
    params.deltaw = state.w_old-state.w;
  end
else
  deltaw_old = params.deltaw;
  params.deltaw = state.w_old-state.w;
  
  params.gamma = params.gamma + (params.deltaw'*deltaw_old) / (deltaw_old'*deltaw_old);
  if params.gamma < 0.5, params.gamma = 0.5; end
end

gamma = params.gamma;
