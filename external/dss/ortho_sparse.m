function [params, result] = ortho_sparse(params, W, w)
% Sparse DSS orthogonalization function
% Assumes that the projection vector has only few non zero items
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
    params.name = 'Sparse orthogonalization';
    params.description = 'Description of this function.';
    return;
end

if nargin>2
  % per component orthogonalization
  if isfield(params,'sparsity')
    sparsity = params.sparsity;
  else
    sparsity = length(w);
  end
  [val,ind] = sort(-abs(w));
  w(ind(sparsity+1:end))=0;
  w = w - W' * W * w;
  result = w / norm(w);
else
  [wdim,vdim] = size(W);
  % symmetric orthogonalization  
  if isfield(params,'sparsity')
    sparsity = params.sparsity;
  else
    sparsity = size(W,2);
  end
  for i = 1 : wdim
    [val,ind] = sort(-abs(W(i,:)));
    W(:,ind(sparsity+1:end))=0;
  end
  W = real(inv(W * W')^(1/2))' * W;
  result = W;
end

%figure(4);
%clf;
%plot(result);
%drawnow;
