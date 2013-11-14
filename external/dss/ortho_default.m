function [params, result] = ortho_default(params, W, w)
% Default DSS orthogonalization function
%   W = orthof(W)     for symmetric dss
%   w = orthof(W, w)  for deflation dss
%     W Matrix with projection vectors as rows. For deflation
%       algorithm only previously calculated projections are given.
%     w Currently iterated projection. Only for deflation algorithm.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
    params.name = 'Default orthogonalization';
    params.description = 'Description of this function.';
    return;
end

if nargin>2
  % per component orthogonalization
  w = w - W' * W * w;
  result = w / norm(w);
else
  % symmetric orthogonalization  
  W = real(inv(W * W')^(1/2))' * W;
  result = W;
end
