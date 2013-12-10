function [A, W] = dss_fastica(Y, wM, dwM, params)

% DSS with FastICA-like interface
%   [A, W] = dss_fastica(Y, wM, dwM, params)
%     Y       Whitened data
%     wM      Whitening matrix
%     dwM     De-whitening matrix
%     params  DSS algorithm parameters (optional)
%     A       Mixing matrix (de-whitened) Y = A * S
%     W       Un-mixing matrix (de-whitened) S = W * Y
%
% Example:
%   [E, D] = pcamat(X);
%   [Y, wM, dwM] = whitenv(X, E, D);
%   [A, W] = dss_fastica(Y, wM, dwM, params);

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if (nargin<4); params=[]; end

% Skip DSS whitening by providing already whitened data
state.Y = Y;
[state, A, W] = denss(state, params);
W = state.W * wM;
A = dwM * state.W';
