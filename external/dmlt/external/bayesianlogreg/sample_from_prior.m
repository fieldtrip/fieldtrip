function S = sample_from_prior(mu, Q, M)
% take M unconditional samples from a GMRF with mean mu and precision matrix Q (algorithm 2.4, Rue)
% samples are generated as column vectors

 % fprintf('sampling from prior\n');

  L = chol(Q,'lower');
  
  if nargin < 2, M = 1; end

  n = length(mu);
  
  v = zeros(n,M);      
  z = randn(n,M);

  % vectorized back substitution
  v(n,:) =  z(n,:)./L(n,n);
  for i=(n-1):-1:1
    v(i,:) = (z(i,:) - (L((i+1):n,i)' * v((i+1):n,:))) ./ L(i,i);
  end

  S = repmat(mu,[1 M]) + v;
