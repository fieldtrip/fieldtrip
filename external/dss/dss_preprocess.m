function state = dss_preprocess(X, state)

% Preprocess data

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if isfield(state, 'Y')
  dss_message(state, 1, 'WARNING: Overwriting existing preprocessed data');
end

if ~isfield(state, 'wdim') state.wdim = size(state.X, 1); end

dss_message(state, 1, 'Preprocessing data\n');
[state.preprocf.params,state.X,state.Xmeans state.Y, state.V, state.dV] = ...
    feval(state.preprocf.h, state.preprocf.params, state.X, state.wdim);

state.wdim = size(state.Y, 1);
if ~isfield(state, 'sdim'); state.sdim = state.wdim; end

dss_message(state,1,sprintf('Preprocessed data dimension %i, extracting  %i components.\n', state.wdim, state.sdim));

