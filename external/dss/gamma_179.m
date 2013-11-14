function [params, gamma] = gamma_179(params, state)
% 179-rule DSS gamma function
%   [params, gamma] = gamma_179(params, state)
%     params  Gamma function parameters
%     state   DSS algorithm state
%     gamma   Result gamma coefficient

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
    params.name = '179-rule';
    params.adaptive=1;
    params.approach = {'pca','defl','symm'};
    return;
end

if state.iteration <= 2
  params.gamma = 1;
  if state.iteration == 2
    params.deltaw = state.w_old-state.w;
  end
elseif params.gamma ~= 1
  deltaw_old = params.deltaw;
  params.deltaw = state.w_old-state.w;

  limit = cos(90/180*pi);
  if params.deltaw'*deltaw_old/norm(params.deltaw)/norm(deltaw_old) <= limit
    fprintf('[gamma = 0.5]');
    params.gamma = 0.5;  
  end
  
end

gamma = params.gamma;

