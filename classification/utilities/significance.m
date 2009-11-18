function [reject,pvalue,level] = significance(cpost1, design1, cpost2, design2, varargin)
% SIGNIFICANCE makes a significance test whether two algorithms perform
% the same or differently
%
% [reject,pvalue,level] = significance(cpost1, design1, cpost2, design2, varargin)
%
% where typically design1 == design2 
%
% reject = reject the null hypothesis? I.e., is there a difference?
% pvalue = used pvalue
% level = used (Bonferroni corrected) significance level
%
% 'level' = 0.05; signifance level for rejecting the null hypothesis
% 'test' = 'mcnemar'; type of significance test
% 'twosided' = true; one versus two-sided test
% 'bonferroni' = 1; lists number of tests; 
%   if larger than 1 we do bonferroni correction for multiple tests
% 'idx'; gives the column index of the design matrix that represents
%    example indices. Used when design1 and design2 are not the same
%
% if inputs are cell arrays then all experiments are
% concatenated (e.g., folds)
% 
% Copyright (C) 2008, Marcel van Gerven
% F.C. Donders Centre for Cognitive Neuroscience, Nijmegen, NL
%
    
    % initialization   

    cfg = varargin2struct(varargin);   

    if ~isfield(cfg,'test'), cfg.test = 'mcnemar'; end
    if ~isfield(cfg,'bonferroni'), cfg.bonferroni = 1; end
    if ~isfield(cfg,'level'), level = 0.05; else level = cfg.level; end
    if ~isfield(cfg,'twosided'), cfg.twosided = true; end    
    
    % change significance level
    if cfg.bonferroni > 1, level = 1 - power(1 - level,1 / cfg.bonferroni); end

    % process cell arrays

    if iscell(cpost1)
        post1 = [];
        for c=1:length(cpost1)
            post1 = [post1; cpost1{c}];
        end
        cpost1 = post1;
    end

    if iscell(cpost2)
        post2 = [];
        for c=1:length(cpost2)
            post2 = [post2; cpost2{c}];
        end
        cpost2 = post2;
    end

    % true classes
    
    if iscell(design1)
        tcls = [];
        for c=1:length(design1)
            tcls = [tcls; design1{c}];
        end
        ctcls1 = tcls;
    else
        ctcls1 = design1;
    end

    if iscell(design2)
        tcls = [];
        for c=1:length(design2)
            tcls = [tcls; design2{c}];
        end
        ctcls2 = tcls;
    else
        ctcls2 = design2;
    end

    % optionally shuffle rows of the design matrices
    if isfield(cfg,'idx')
       
        [a,b] = sort(ctcls1(:,cfg.idx));        
        ctcls1 = ctcls1(b,:);
        cpost1 = cpost1(b,:);
        
        [a,b] = sort(ctcls2(:,cfg.idx));        
        ctcls2 = ctcls2(b,:);
        cpost2 = cpost2(b,:);
       
        % now both designs should be comparable
    else
       
        % check if designs are the same if index is left unspecified
        if any(ctcls1(:,1) ~= ctcls2(:,1))
            error('designs do not match; please specify indices');
        end
    end

    % execute a one-sided test
    switch cfg.test

        case 'mcnemar'
            pvalue = mcnemar(cpost1,cpost2,ctcls1);

    end

    % make two-sided if necessary
    if cfg.twosided, level = level/2; end

    % make decision if we can reject the null hypothesis
    reject = (pvalue < level);

end

function pvalue = mcnemar(cpost1, cpost2, ctcls2)
% MCNEMAR tests whether the posterior probabilities cpost1 and cpost2
% are significantly different.
%
% Note: this is an approximation to the computationally heavy binomial
% test. Later, this file should become a general significance testing
% function that also accepts cell arrays for cpost1 and cpost2

    % compute class labels
    [temp,pcls1] = max(cpost1,[],2);
    [temp,pcls2] = max(cpost2,[],2);

    % compute winners for each algorithm

    p1 = (pcls1 == ctcls2(:,1));
    p2 = (pcls2 == ctcls2(:,1));

    s = sum(p1 & ~p2); % successes w.r.t. 1
    f = sum(~p1 & p2); % failures w.r.t. 1

    % compute one sided p-values

    if (s+f)
      mcnemarstat = (abs(s - f) - 1).^2 / (s+f);
    else
      mcnemarstat = 0;
    end

    pvalue = 1 - chi2cdf(mcnemarstat,1);

end
