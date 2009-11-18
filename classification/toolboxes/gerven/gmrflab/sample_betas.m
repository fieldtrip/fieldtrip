function [S,u,v] = sample_betas(Q, M)
% take M unconditional samples from Bayesian logistic regression model
% where the auxiliary variables u and v have prior precision matrix Q

  fprintf('sampling betas from auxiliary variables\n');
 
  if nargin < 2, M = 1; end
  n = size(Q,1);
 
  % get samples for auxiliary variables
  u = sample_from_prior(zeros(n,1),Q,M);
  v = sample_from_prior(zeros(n,1),Q,M);

  % get samples for betas
  S = normrnd(zeros(n,M),sqrt(u.^2 + v.^2));
