function pvalue = ft_mv_mcnemar(design,post1,post2)
% FT_MV_MCNEMAR performs a mcnemar test which is an approximate binomial
% test. It compares if post1 is significantly better than post2 given the real outcome (design). 
%
% ft_mv_mcnemar(design,post1,post2)
%
% RETURNS
% pvalue      estimated pvalue
%
% Ref:
% Salzberg, On Comparing Classifiers: Pitfalls to Avoid and a Recommended Approach
%
% Copyright (c) 2010, Marcel van Gerven

  % compute p-value using one-sided (approximate) mcnemar test
  
  % compute class labels
  [tmp,pcls1] = max(post1,[],2);
  [tmp,pcls2] = max(post2,[],2);
  
  % compute winners for each algorithm
  p1 = (pcls1 == design);
  p2 = (pcls2 == design);
  
  s = sum(p1 & ~p2); % successes w.r.t. 1
  f = sum(~p1 & p2); % failures w.r.t. 1
  
  % compute one-tailed mcnemar test
  if (s+f) > 0 && s > f
    mcnemarstat = (abs(s - f) - 1).^2 / (s+f);
    pvalue = 1 - chi2cdf(mcnemarstat,1); % two-tailed
    pvalue = pvalue / 2; % one-tailed
  else
    pvalue = 1;
  end
  
  
end