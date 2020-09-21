function [params,X,means,wX,wM,dwM] = pre_sphere(params, X, dim)
% Spheres the data with the eigenvalue decomposition of the covariance matrix
%   [params,X,means,wX,wM,dwM] = pre_sphere(params, X, dim)
%     X    Data to be sphered
%     dim  Requested result data dimension
%     wX   Sphered data
%     wM   Sphering matrix
%     dwM  De-sphering matrix

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
    params.name = 'PCA sphering';
    params.description = 'Spheres the data with the eigenvalue decomposition of the covariance matrix.'; 
    return;
end

if nargin<3; dim=0; end

[xdim,tdim] = size(X);
% removing the mean
means = mean(X,2);
X = X - repmat(means,1,tdim);

[wM, dwM] = dss_sphere(X, dim);

wX = wM * X;

