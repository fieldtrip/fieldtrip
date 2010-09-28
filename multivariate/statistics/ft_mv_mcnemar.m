function [pvalue,reject,level] = ft_mv_mcnemar(design,post1,post2,varargin)
% FT_MV_MCNEMAR performs a mcnemar test which is an approximate binomial
% test. It compares if two predicted posteriors (post1 and post2) are significantly different 
% given the real outcome (design). 
%
% OPTIONS
% bonferroni  Bonferroni correction, i.e., nr of comparisons (1)
% level       significance level (0.05)
% twosided    twosided test (true)
%
% RETURNS
% reject      whether or not the null hypothesis that both outcomes are the same is rejected
% pvalue      estimated pvalue
% level       used threshold after Bonferroni correction
%
% Ref:
% Salzberg, On Comparing Classifiers: Pitfalls to Avoid and a Recommended Approach
%
% Copyright (c) 2010, Marcel van Gerven

  cfg = struct(varargin{:});
  if ~isfield(cfg,'bonferroni'), cfg.bonferroni = 1; end
  if ~isfield(cfg,'level'), level = 0.05; else level = cfg.level; end
  if ~isfield(cfg,'twosided'), cfg.twosided = true; end
  
  % bonferroni correction
  if cfg.bonferroni > 1, level = 1 - power(1 - level,1 / cfg.bonferroni); end
    
  % compute p-value using one-sided (approximate) mcnemar test
  
  % compute class labels
  [tmp,pcls1] = max(post1,[],2);
  [tmp,pcls2] = max(post2,[],2);
  
  % compute winners for each algorithm
  
  p1 = (pcls1 == design);
  p2 = (pcls2 == design);
  
  s = sum(p1 & ~p2); % successes w.r.t. 1
  f = sum(~p1 & p2); % failures w.r.t. 1
  
  % compute one sided p-values
  if (s+f)
    mcnemarstat = (abs(s - f) - 1).^2 / (s+f);
  else
    mcnemarstat = 0;
  end
  pvalue = 1 - chi2cdf(mcnemarstat,1);
  
  % make two-sided if necessary
  if cfg.twosided, level = level/2; end
  
  % make decision if we can reject the null hypothesis
  reject = (pvalue < level);
  
end