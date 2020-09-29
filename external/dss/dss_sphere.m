function [V,dV] = dss_sphere(X, dim, symmetric, params)
%   [V,dV] = dss_sphere(X)
%   [V,dV] = dss_sphere(X, dim)
%   [V,dV] = dss_sphere(X, dim, symmetric)
%     X          Data to be sphered
%     dim        Requested result data dimension, 0 for no reduction
%     symmetric  Boolean, perform symmetric sphering
%     V          Sphering matrix
%     dV         De-sphering matrix

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: dss_sphere.m,v 1.2 2005/11/30 08:29:40 jaakkos Exp $

if nargin<4
  params = [];
end
if nargin<3 || isempty(symmetric)
  symmetric = false;
end
if ~isfield(params, 'absthreshold')
  params.absthreshold=true;
end
if params.absthreshold && ~isfield(params, 'thresholdEig')
  % Treshold for eigenvalues that are treated zero
  params.thresholdEig=1e-12;
elseif ~params.absthreshold && ~isfield(params, 'thresholdEig')
  params.thresholdEig=0.99; % keep the components that explain up to 99% of variance
end
if ~isfield(params, 'thresholdSphered')
  % Treshold for maximum covariance for treating data already sphered
  params.tresholdSphered=1e-12;
end

if nargin>=2 && dim~=0
  if iscell(X)
    dim = min(dim, size(X{1}, 1));
  else
    dim = min(dim, size(X, 1));
  end 
elseif iscell(X)
  dim = size(X{1}, 1);
else
  dim = size(X, 1);
end

fprintf('Sphering the data\n');
% data covariance
if iscell(X)
  C = cellcov(X, 2, 1);
  I = diag(ones(size(X{1},1),1));
else
  C = cov(X', 1);
  I = diag(ones(size(X,1),1));
end

% check if data is already sphered
maxCov = max(max(C-I));
if maxCov<params.tresholdSphered
    fprintf('SPHERING: Maximum data covariance %d below treshold (%d), assume pre-sphered data.\n', maxCov, params.tresholdSphered);
    V = diag(ones(dim,1));
    dV = V;
else

    % eigenvectors and eigenvalues
    [E,D] = eig(C);
    D = diag(D);

    % sort eigenvalues in descending order
    [D,order] = sort(-D);
    D=-D;
    E=E(:,order);

    % reduce dimension based on matrix rank
    if params.absthreshold
      values = sum(D>=params.thresholdEig);
      if values<length(D)
        fprintf('SPHERING: %d values (%d negative) below eigenvalue treshold (%d).\n', length(D)-values, sum(D<0), params.thresholdEig);
      end
    else
      values = sum(cumsum(D)./sum(D)<=params.thresholdEig);
      if values<length(D)
        fprintf('SPHERING: %d values (%d negative) below cumulative eigenvalue treshold (%d).\n', length(D)-values, sum(D<0), params.thresholdEig);
      end
    end 
    % reduce dimension
    dim = min(dim, values);
    if iscell(X) && dim < size(X{1},1)
      fprintf('SPHERING: Reducing data dimension: %d -> %d.\n', size(X{1},1), dim);
    end
    if ~iscell(X) && dim < size(X,1)
      fprintf('SPHERING: Reducing data dimension: %d -> %d.\n', size(X,1), dim);
    end
    E = E(:,1:dim);
    D = D(1:dim);

    % calculate sphering matrices
    D = sqrt(D);
    V = diag(1./D) * E' ;
    dV = E * diag(D);
    if symmetric
        % symmetric sphering
        V = E * V;
        dV = dV *E';
    end

end
