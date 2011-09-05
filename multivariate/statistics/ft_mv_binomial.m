function pvalue = ft_mv_binomial(design,post1,post2,varargin)
% FT_MV_BINOMIAL performs a binomial test. It compares if post1 is significantly better than post2 given the real outcome (design). 
%
% pvalue = ft_mv_binomial(design,post1,post2,varargin)
%
% RETURNS
% pvalue      estimated pvalue
%
% Ref:
% Salzberg, On Comparing Classifiers: Pitfalls to Avoid and a Recommended Approach
%
% Copyright (c) 2010, Marcel van Gerven

  % compute class labels
  [tmp,pcls1] = max(post1,[],2);
  [tmp,pcls2] = max(post2,[],2);
  
  % compute winners for each algorithm
  p1 = (pcls1 == design);
  p2 = (pcls2 == design);
  
  s = sum(p1 & ~p2); % successes w.r.t. 1
  f = sum(~p1 & p2); % failures w.r.t. 1
  
  % compute one-tailed binomial test
  if (s+f) > 0 && s > f % if there's a difference and it's in favor of p1 over p2
    pvalue = 1 - binocdf(s-1,s+f,0.5);
  else
    pvalue = 1;
  end
  
end