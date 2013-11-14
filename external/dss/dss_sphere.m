function [V,dV] = dss_sphere(X, dim, symmetric)
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
% $Id$

% Treshold for eigenvalues that are treated zero
tresholdEig=1e-9;
% Treshold for maximum covariance for treating data already sphered
tresholdSphered=1e-9;

if nargin>=2 & dim~=0; dim = min(dim, size(X, 1));
else; dim = size(X, 1);
end

fprintf('Sphering the data\n');
% data covariance
C = cov(X', 1);

% check if data is already sphered
I = diag(ones(size(X,1),1));
maxCov = max(max(C-I));
if maxCov<tresholdSphered
    fprintf('SPHERING: Maximum data covariance %d below treshold (%d), assume pre-sphered data.\n', maxCov, tresholdSphered);
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
    values = sum(D>=tresholdEig);
    if values<length(D)
        fprintf('SPHERING: %d values (%d negative) below eigenvalue treshold (%d).\n', length(D)-values, sum(D<0), tresholdEig);
    end

    % reduce dimension
    dim = min(dim, values);
    if (dim < size(X,1))
      fprintf('SPHERING: Reducing data dimension: %d -> %d.\n', size(X,1), dim);
    end
    E = E(:,1:dim);
    D = D(1:dim);

    % calculate sphering matrices
    D = sqrt(D);
    V = diag(1./D) * E' ;
    dV = E * diag(D);
    if nargin>=3 & symmetric
        % symmetric sphering
        V = E * V;
        dV = dV *E';
    end

end
