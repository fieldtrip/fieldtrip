function state = dss_core_pca(state)
% PCA DSS algorithm
%   state = dss_core_pca(state)
%     Calculates linear denoising source separation with PCA.
%     Performs calculation defined by state structure and
%     returns the result as state structure.
%     Used CPU time is recorded in the state.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: dss_core_pca.m,v 1.13 2005/12/02 12:23:18 jaakkos Exp $

dss_message(state,2,'Extracting components in PCA DSS\n');

% -- Initialize local variables
start_time = cputime;
interrupt_iteration = 0;
interrupt_component = 0;
% Dimension of the whitened data
if iscell(state.Y)
  wdim = size(state.Y{1},1);
else
  wdim = size(state.Y, 1);
end
sdim = state.sdim;

% -- Denoising
[state.denf.params, Y] = feval(state.denf.h, state.denf.params, state.Y, state);

% -- PCA
if iscell(Y)
  [E,D] = eig(cellcov(Y,2,1));
else
  [E,D] = eig(Y * Y' / size(Y,2));
end
D = diag(D);

% sort eigenvalues in descending order
[D,order] = sort(-D);
D=-D;
E=E(:,order(1:min(sdim,size(E,2))));

state.D = D;
state.W = E';
state.S = state.W * state.Y;

% -- record the used cpu time
if ~isfield(state, 'cputime'); state.cputime = 0; end
state.cputime = state.cputime + cputime - start_time;
