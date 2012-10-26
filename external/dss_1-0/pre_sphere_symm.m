function [params,X,means,wX,wM,dwM] = pre_sphere_symm(params, X, dim)
% Spheres the data with covarianve matrix eigenvalue decomposition.
%   [params,X,means,wX,wM,dwM] = pre_sphere_symm(params, X, dim)
%     X    Data to be sphered
%     dim  Requested result data dimension
%     wX   Sphered data
%     wM   Sphering matrix
%     dwM  De-sphering matrix

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: pre_sphere_symm.m,v 1.7 2005/11/30 08:29:40 jaakkos Exp $

if nargin<2
    params.name = 'Symmetric sphering';
    params.description = 'Spheres the data with covarianve matrix eigenvalue decomposition.';
    return;
end

if nargin<3; dim=0; end

[xdim,tdim] = size(X);
% removing the mean
means = mean(X,2);
X = X - repmat(means,1,tdim);

[wM, dwM] = dss_sphere(X, dim, 1);

wX = wM * X;
