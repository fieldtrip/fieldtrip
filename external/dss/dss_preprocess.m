function state = dss_preprocess(X, state)

% Preprocess data

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: dss_preprocess.m,v 1.6 2005/11/30 08:29:40 jaakkos Exp $

if isfield(state, 'Y')
  dss_message(state, 1, 'WARNING: Overwriting existing preprocessed data');
end

if ~isfield(state, 'wdim')
  if iscell(state.X)
    state.wdim = size(state.X{1},1);
  else
    state.wdim = size(state.X, 1);
  end
end

dss_message(state, 1, 'Preprocessing data\n');
[state.preprocf.params,state.X,state.Xmeans state.Y, state.V, state.dV] = ...
    feval(state.preprocf.h, state.preprocf.params, state.X, state.wdim);

if iscell(state.Y)
  state.wdim = size(state.Y{1},1);
else
  state.wdim = size(state.Y, 1);
end

% broadcast the block indx (if existing) to the state
if isfield(state.preprocf.params, 'indx')
  state.indx = state.preprocf.params.indx;
end

if ~isfield(state, 'sdim'); state.sdim = state.wdim; end

dss_message(state,1,sprintf('Preprocessed data dimension %i, extracting  %i components.\n', state.wdim, state.sdim));

