function state = dss_core_mim(state)
% MIM DSS algorithm
%   state = dss_core_mim(state)
%     Calculates linear denoising source separation with MIM.
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
  C = cellcov(Y,2,0);
  %[E,D] = eig(cellcov(Y,2,0));
else
  C = (Y*Y')/size(Y,2);
  %[E,D] = eig(Y * Y' / size(Y,2));
end

i1 = state.indx==1;
i2 = state.indx==2;

% create whitening matrix T
T = C;
T(i1, i2) = 0;
T(i2, i1) = 0;
T = inv(sqrtm(real(T)));

E = imag(T(i1, i1) * C(i1, i2) *T(i2, i2)');

[e1, d1] = eig(E*E');
[e2, d2] = eig(E'*E);

d1 = diag(d1);
d2 = diag(d2);
n  = min(sum(i1), sum(i2));

% sort eigenvalues in descending order
[d1,order] = sort(-d1);
d1 = -d1(1:n);
e1 = e1(:,order(1:n));

[d2,order] = sort(-d2);
d2 = -d2(1:n);
e2 = e2(:,order(1:n));


state.D = [d1;d2];
state.sdim = numel(state.D);
state.W = blkdiag(e1,e2)'/T;
state.S = state.W * state.Y;

% -- record the used cpu time
if ~isfield(state, 'cputime'); state.cputime = 0; end
state.cputime = state.cputime + cputime - start_time;
