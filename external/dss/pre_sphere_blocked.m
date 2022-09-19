function [params,X,means,wX,wM,dwM] = pre_sphere_blocked(params, X, dim)
% Spheres the data with the eigenvalue decomposition of the covariance matrix
%   [params,X,means,wX,wM,dwM] = pre_sphere_blocked(params, X, dim)
%     X    Data to be sphered
%     dim  Requested result data dimension
%     wX   Sphered data
%     wM   Sphering matrix
%     dwM  De-sphering matrix

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: pre_sphere.m,v 1.18 2005/11/30 08:20:52 jaakkos Exp $

if nargin<2
    params.name = 'PCA sphering';
    params.description = 'Spheres the data with the eigenvalue decomposition of the covariance matrix.'; 
    return;
end

if nargin<3; dim=0; end

if isfield(params, 'indx') 
  indx = params.indx;
end

[xdim,tdim] = size(X);
% removing the mean
means = mean(X,2);
if iscell(X)
  X = cellvecadd(X, -means);
  if isempty(indx)
    indx = ones(1,size(X{1},1));
  end
  uindx = unique(indx);
  newindx = zeros(1,0);
  for k = 1:numel(uindx)
    ix = indx==uindx(k);
    [wmat, dwmat] = dss_sphere(cellrowselect(X,ix), sum(ix), [], params);
    newindx       = cat(2, newindx, ones(1,size(wmat,1)).*uindx(k));
    iy = newindx==uindx(k);
    
    wM(iy,ix)  = wmat;
    dwM(ix,iy) = dwmat;
  end
  params.indx = newindx;
  
else
  X = X - repmat(means,1,tdim);
  if isempty(indx)
    indx = ones(1,size(X,1));
  end
  uindx = unique(indx);
  for k = 1:numel(uindx)
    ix = indx==uindx(k);
    
    [wM(ix,ix), dwM(ix,ix)] = dss_sphere(X(ix,:), sum(ix), [], params);
  end
  
end

wX = wM * X;

