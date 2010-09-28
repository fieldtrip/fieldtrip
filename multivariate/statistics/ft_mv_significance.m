function res = ft_mv_significance(design,post,sigtest)
% FT_MV_SIGNIFICANCE signficance tests for a real outcome (design) and predicted
% outcome (post) and returns a p-value
%
% Copyright (c) 2010, Marcel van Gerven

  if isempty(sigtest)
    error('significance test undefined');
  end
  
  if ~iscell(design)
    design = {design};
    post = {post};
  end
  
  res = cell(length(design),1);
  for c=1:length(design)
    
    switch lower(sigtest)
      
      case 'mcnemar'
        % one-sided approximate binomial test which compares outcome with that of a majority classifier
        
        % create majority classifier results
        nclasses = size(post{c},2);
        priors = zeros(1,nclasses);
        for k=1:nclasses
          priors(k) = priors(k) + sum(design{c}(:,1)==k);
        end
        [tmp,clss] = max(priors);
        rndpost = zeros(size(post{c}));
        rndpost(:,clss) = 1;
        
        res{c} = ft_mv_mcnemar(design{c},post{c},rndpost,'twosided',false);
        
    end
    
  end
  
  if length(res)==1, res = res{1}; end
  
end

