function [w,rho] = bsscca(X, delay)

% BSSCCA computes the mixing matrix based on the canonical correlation between a signal and its lagged-one copy. It implements the algorithm described in [1]
%
% DeClercq et al 2006, IEEE Biomed Eng 2583.

if nargin<2,
  delay = 1;
end

[n,m] = size(X);

% get the means
%m   = ones(1,m-1);
%mX  = mean(X(:,2:end),2);   % lag zero
%mX2 = mean(X(:,1:end-1),2); % lag one

% use Borga's (2001) formulation from 'a unified approach to PCA, PLS, MLR
% and CCA'
A = zeros(2*n);
B = zeros(2*n);

XY = X(:,delay+(1:m))*X(:,1:(m-delay))';
%XY = (X(:,2:end)-mX*m)*(X(:,1:end-1)-mX2*m)';

A(1:n,(n+1):end) = XY;
A((n+1):end,1:n) = XY';
%B(1:n,1:n)       = (X(:,2:end)-mX*m)*(X(:,2:end)-mX*m)';
%B((n+1):end,(n+1):end) = (X(:,1:end-1)-mX2*m)*(X(:,1:end-1)-mX2*m)';
B(1:n,1:n)       = X(:,delay+(1:m))*X(:,delay+(1:m))';
B((n+1):end,(n+1):end) = X(:,1:(m-delay))*X(:,1:(m-delay))';

[w,rho]  = eig(B\A);
[~,ix]   = sort(diag(rho.^2),'descend');
w        = w(1:n,ix(2:2:end))';
rho      = rho(ix(2:2:end),ix(2:2:end)).^2;

% normalise to unit norm
for k= 1:size(w,1)
  w(k,:) = w(k,:)./norm(w(k,:));
end
