function [params, result] = ortho_quasi(params, W, w)
% Quasi-orthogonalization
%   W = orthof(params, W)     for symmetric dss
%   w = orthof(params, W, w)  for deflation dss
%     params.alpha  For deflation...
%     W Matrix with projection vectors as rows. For deflation
%       algorithm only previously calculated projections are given.
%     w Currently iterated projection. Only for deflation algorithm.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
    params.name = 'Quasi-orthogonalization';
    params.description = 'Description of this function.';
    params.param = {'alpha'};
    params.param_type = {'scalar'};
    params.param_value = {0.1};
    params.param_desc = {'For deflation...'};
    return;
end

if ~isfield(params, 'alpha')
  % TODO: set alpha based on data dimension
  params.alpha = 0.1;
end

if nargin>2
  % per component quasi-orthogonalization
  w = w - params.alpha * W' * W * w;
  result = w / norm(w);
else
  % symmetric quasi-orthogonalization  
  W = 3/2*W - W'*W*W'/2;
  result = W.*(sum(W.^2,2).^(-1/2)*ones(1,size(W,2)));
end
